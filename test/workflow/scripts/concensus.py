import argparse
import sys
from collections import defaultdict
from contextlib import ExitStack

import cigar
import polars as pl
import pysam
from intervaltree import Interval, IntervalTree

DEF_PAF_COLS = (
    "qname",
    "qlen",
    "qst",
    "qend",
    "ort",
    "tname",
    "tlen",
    "tst",
    "tend",
    "matches",
    "aln_blk_len",
    "q",
    "NM",
    "ms",
    "AS",
    "nn",
    "tp",
    "cg",
)

CTG_COORD_COLS = [
    "qname_ctg_start",
    "tname_ctg_start",
    "qname_ctg_end",
    "tname_ctg_end",
]


def create_intervaltree(file: str) -> defaultdict[str, IntervalTree]:
    intervals = defaultdict(IntervalTree)
    df_query_misasm = pl.read_csv(
        file,
        has_header=False,
        separator="\t",
        new_columns=["name", "start", "stop", "misassembly"],
    )
    for row in df_query_misasm.iter_rows():
        intervals[row[0]].add(Interval(row[1] - 1, row[2] - 1, row[3]))

    return intervals


def main():
    ap = argparse.ArgumentParser(
        description="Generate a concensus sequence from an alignment between two contigs."
    )
    ap.add_argument(
        "-i",
        "--input_paf",
        type=str,
        required=True,
        help="Alignment PAF broken via rustybam break-paf.",
    )
    ap.add_argument(
        "-r",
        "--input_ref_fa",
        type=str,
        required=True,
        help="Reference fasta sequence.",
    )
    ap.add_argument(
        "-q", "--input_query_fa", type=str, required=True, help="Query fasta sequence."
    )
    ap.add_argument(
        "-o",
        "--output_fa",
        default=sys.stdout,
        type=argparse.FileType("wt"),
        help="Output fasta. Can be merged with the --merge flag.",
    )
    ap.add_argument(
        "--input_ref_misasm_bed",
        type=str,
        required=True,
        help="Reference misassembly bed.",
    )
    ap.add_argument(
        "--input_query_misasm_bed",
        type=str,
        required=True,
        help="Query misassembly bed.",
    )
    ap.add_argument(
        "-b",
        "--output_bed",
        default=None,
        type=str,
        help="Output bedfile.",
    )
    ap.add_argument("--merge", action="store_true", help="Merge output_fa sequences.")
    ap.add_argument("--verbose", action="store_true", help="Verbose logging output.")

    args = ap.parse_args()

    with ExitStack() as es:
        ref_fa = es.enter_context(pysam.FastaFile(args.input_ref_fa))
        query_fa = es.enter_context(pysam.FastaFile(args.input_query_fa))

        df_ref_lens = pl.DataFrame(
            [(ctg, length) for ctg, length in zip(ref_fa.references, ref_fa.lengths)],
            orient="row",
            schema={"base_tname": pl.String, "length": pl.Int64},
        )
        df_query_lens = pl.DataFrame(
            [
                (ctg, length)
                for ctg, length in zip(query_fa.references, query_fa.lengths)
            ],
            orient="row",
            schema={"base_qname": pl.String, "length": pl.Int64},
        )
        df_paf = (
            pl.read_csv(
                args.input_paf,
                separator="\t",
                has_header=False,
            )
            .rename({f"column_{i}": c for i, c in enumerate(DEF_PAF_COLS, 1)})
            .with_columns(
                base_qname=pl.col("qname").str.split_exact(":", 1),
                base_tname=pl.col("tname").str.split_exact(":", 1),
            )
            .filter((pl.col("tp") == "tp:A:P") & (pl.col("q") == 60))
            .unnest("base_qname")
            .rename({"field_0": "base_qname", "field_1": "qname_coords"})
            .unnest("base_tname")
            .rename({"field_0": "base_tname", "field_1": "tname_coords"})
            .with_columns(
                pl.col("qname_coords").str.split_exact("-", 1),
                pl.col("tname_coords").str.split_exact("-", 1),
            )
            .unnest("qname_coords")
            .rename({"field_0": "qname_ctg_start", "field_1": "qname_ctg_end"})
            .unnest("tname_coords")
            .rename({"field_0": "tname_ctg_start", "field_1": "tname_ctg_end"})
            .cast({col: pl.Int64 for col in CTG_COORD_COLS})
            .with_columns(
                qname=pl.col("base_qname")
                + ":"
                + (pl.col("qname_ctg_start") - 1).cast(pl.String)
                + "-"
                + pl.col("qname_ctg_end").cast(pl.String),
                tname=pl.col("base_tname")
                + ":"
                + (pl.col("tname_ctg_start") - 1).cast(pl.String)
                + "-"
                + pl.col("tname_ctg_end").cast(pl.String),
            )
            .with_columns(
                pl.col("tst") + pl.col("tname_ctg_start"),
                pl.col("tend") + pl.col("tname_ctg_start"),
                pl.col("qst") + pl.col("qname_ctg_start"),
                pl.col("qend") + pl.col("qname_ctg_start"),
            )
            .join(
                df_query_lens,
                on=["base_qname"],
            )
            .rename({"length": "query_ctg_len"})
            .join(
                df_ref_lens,
                on=["base_tname"],
            )
            .rename({"length": "target_ctg_len"})
            .sort(by="qname")
            .with_row_index()
        )

        target_misasm_intervals = create_intervaltree(args.input_ref_misasm_bed)
        query_misasm_intervals = create_intervaltree(args.input_query_misasm_bed)

        new_ctgs = []
        primary_alns = df_paf.with_columns(
            pl.col("cg").str.replace("cg:Z:", "", literal=True)
        ).iter_rows(named=True)
        for aln in primary_alns:
            cg = cigar.Cigar(aln["cg"])
            target_bp_added = 0
            query_bp_added = 0
            for bp, op in cg.items():
                adj_target_start = aln["tst"] + target_bp_added
                adj_query_start = aln["qst"] + query_bp_added
                target_interval = Interval(adj_target_start, adj_target_start + bp)
                if op == "=":
                    new_ctgs.append(
                        (
                            aln["base_tname"],
                            "tname",
                            adj_target_start,
                            adj_target_start + bp,
                            aln["target_ctg_len"],
                        )
                    )
                    target_bp_added += bp
                    query_bp_added += bp
                # Deletion in query with respect to reference.
                elif op == "D":
                    # # TODO: Will depend on if false-duplication or gap in sequence.
                    # new_ctgs.append(("", "spacer", 1, 1, 1))

                    overlapping_target_misasms = target_misasm_intervals[
                        aln["tname"]
                    ].overlap(target_interval)
                    overlapping_query_misasms = query_misasm_intervals[
                        aln["qname"]
                    ].overlap(target_interval)

                    if args.verbose:
                        print("# Deletion in query relative to target:", target_interval, file=sys.stderr)
                        print(
                            "Target name:",
                            aln["tname"],
                            overlapping_target_misasms,
                            file=sys.stderr,
                        )
                        print(
                            "Query name:",
                            aln["qname"],
                            overlapping_query_misasms,
                            file=sys.stderr,
                        )

                    if overlapping_target_misasms and overlapping_query_misasms:
                        largest_target_misassembly: Interval = max(overlapping_target_misasms, key=lambda x: x.begin - x.end)
                        largest_query_misassembly: Interval = max(overlapping_query_misasms, key=lambda x: x.begin - x.end)
                        match (largest_target_misassembly.data, largest_query_misassembly.data):
                            # IF was just gap, could conclude that target had sequence added.
                            # BUT because MISJOIN, deletion in query explained.
                            case ("GAP", "MISJOIN"):
                                pass
                            case _:
                                pass
                    elif overlapping_target_misasms:
                        largest_misassembly: Interval = max(overlapping_target_misasms, key=lambda x: x.begin - x.end)
                        match largest_misassembly.data:
                            # If misjoin or gap in query, sequence in target is added. Remove it from target.
                            case "MISJOIN" | "GAP":
                                new_ctgs.append(("", "spacer", 1, 1, 1))
                            case "COLLAPSE":
                                pass
                            case "COLLAPSE_VAR":
                                pass
                            case _:
                                pass
                        pass
                    elif overlapping_query_misasms:
                        largest_misassembly: Interval = max(overlapping_query_misasms, key=lambda x: x.begin - x.end)

                        pass

                    target_bp_added += bp
                # Insertion in query with respect to reference.
                elif op == "I":
                    overlapping_target_misasms = target_misasm_intervals[
                        aln["tname"]
                    ].overlap(target_interval)
                    overlapping_query_misasms = query_misasm_intervals[
                        aln["qname"]
                    ].overlap(target_interval)
                    if args.verbose:
                        print("# Insertion in query relative to target:", target_interval, file=sys.stderr)
                        print(
                            "Target name:",
                            aln["tname"],
                            overlapping_target_misasms,
                            file=sys.stderr,
                        )
                        print(
                            "Query name:",
                            aln["qname"],
                            overlapping_query_misasms,
                            file=sys.stderr,
                        )
                    if overlapping_target_misasms and overlapping_query_misasms:
                        pass
                    elif overlapping_target_misasms:
                        largest_misassembly: Interval = max(overlapping_target_misasms, key=lambda x: x.begin - x.end)
                        match largest_misassembly.data:
                            # If misjoin or gap, sequence in target is missing. Add it from query.
                            case "MISJOIN" | "GAP":
                                new_ctgs.append(
                                    (
                                        aln["base_qname"],
                                        "qname",
                                        adj_query_start,
                                        adj_query_start + bp,
                                        aln["query_ctg_len"],
                                    )
                                )
                            case "COLLAPSE":
                                pass
                            case "COLLAPSE_VAR":
                                pass
                            case _:
                                pass
                        pass
                    elif overlapping_query_misasms:
                        largest_misassembly = max(overlapping_query_misasms, key=lambda x: x.begin - x.end)
                        match largest_misassembly.data:
                            case "MISJOIN" | "GAP":
                                pass
                            case "COLLAPSE":
                                pass
                            case "COLLAPSE_VAR":
                                pass
                            case _:
                                pass

                    query_bp_added += bp
                else:
                    target_bp_added += bp
                    query_bp_added += bp
            # TODO: Gaps - check next row.
            try:
                next_aln = df_paf.row(aln["index"] + 1, named=True)
                new_ctgs.append(
                    (
                        aln["base_tname"],
                        "tname",
                        aln["tend"],
                        next_aln["tst"],
                        aln["target_ctg_len"],
                    )
                )
            except pl.exceptions.OutOfBoundsError:
                pass

        df_new_ctgs = pl.DataFrame(
            new_ctgs,
            orient="row",
            schema={
                "ctg_name": pl.String,
                "name": pl.String,
                "st": pl.Int64,
                "end": pl.Int64,
                "ctg_len": pl.Int64,
            },
        )
        df_new_ctgs = (
            df_new_ctgs
            # Group and aggregate to min max coords
            .with_columns(row_grp_id=pl.col("name").rle_id())
            .group_by(["row_grp_id"])
            .agg(
                pl.col("ctg_name").first(),
                pl.col("name").first(),
                pl.col("st").min(),
                pl.col("end").max(),
                pl.col("ctg_len").first(),
            )
            .sort(by="row_grp_id")
            # Remove group agg spacers.
            .filter(pl.col("name") != "spacer")
            .drop("row_grp_id")
            # Then add boundary coordinates.
            .with_row_index()
            .with_columns(
                st=pl.when(pl.col("index") == pl.col("index").first())
                .then(pl.lit(0))
                .otherwise(pl.col("st")),
                end=pl.when(pl.col("index") == pl.col("index").last())
                .then(pl.col("ctg_len"))
                .otherwise(pl.col("end")),
            )
            .drop("index", "ctg_len")
        )

        if df_new_ctgs.is_empty():
            return

        if args.merge:
            header = df_new_ctgs["ctg_name"][0]
            args.output_fa.write(f">{header}\n")

        for row in df_new_ctgs.iter_rows(named=True):
            fa_file = query_fa if row["name"] == "qname" else ref_fa
            if not args.merge:
                new_header = f'>{row["ctg_name"]}:{row["st"]}-{row["end"]}\n'
                args.output_fa.write(new_header)
                args.output_fa.write(
                    fa_file.fetch(row["ctg_name"], start=row["st"], end=row["end"])
                )
                args.output_fa.write("\n")
            else:
                args.output_fa.write(
                    fa_file.fetch(row["ctg_name"], start=row["st"], end=row["end"])
                )

        if args.output_bed:
            df_new_ctgs.write_csv(args.output_bed, separator="\t", include_header=False)


if __name__ == "__main__":

    raise SystemExit(main())
