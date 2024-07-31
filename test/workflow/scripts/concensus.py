import argparse
import sys
from contextlib import ExitStack

import polars as pl
import pysam

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
    "id",
    "cg",
)

CTG_COORD_COLS = [
    "qname_ctg_start",
    "tname_ctg_start",
    "qname_ctg_end",
    "tname_ctg_end",
]


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
    ap.add_argument("--merge", action="store_true", help="Merge output_fa sequences.")

    args = ap.parse_args()

    with ExitStack() as es:
        ref_fa = es.enter_context(pysam.FastaFile(args.input_ref_fa))
        query_fa = es.enter_context(pysam.FastaFile(args.input_query_fa))

        df_ref_lens = pl.DataFrame(
            [(ctg, length) for ctg, length in zip(ref_fa.references, ref_fa.lengths)],
            schema={"base_tname": pl.String, "length": pl.Int64},
        )
        df_query_lens = pl.DataFrame(
            [
                (ctg, length)
                for ctg, length in zip(query_fa.references, query_fa.lengths)
            ],
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
        # pl.Config(tbl_cols=21, tbl_width_chars=1000, fmt_str_lengths=1000)

        new_rows = []
        for row in df_paf.iter_rows(named=True):
            new_rows.append((row["base_qname"], "qname", row["qst"], row["qend"], row["query_ctg_len"]))
            try:
                next_row = df_paf.row(row["index"] + 1, named=True)
            except pl.OutOfBoundsError:
                continue
            if (
                row["qend"] == next_row["qst"]
                and row["q"] == next_row["q"]
                and row["ort"] == next_row["ort"]
            ):
                # Then adjust orientation.
                if row["ort"] == "+":
                    new_st, new_end = row["tend"], next_row["tst"]
                else:
                    new_st, new_end = (
                        row["target_ctg_len"] - row["tend"],
                        row["target_ctg_len"] - next_row["tst"],
                    )
                new_rows.append((row["base_tname"], "tname", new_st, new_end, row["target_ctg_len"]))

        df_new_rows = (
            pl.DataFrame(
                new_rows,
                schema={
                    "ctg_name": pl.String,
                    "name": pl.String,
                    "st": pl.Int64,
                    "end": pl.Int64,
                    "ctg_len": pl.Int64,
                },
            )
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
            .drop("index")
        )

        if df_new_rows.is_empty():
            return

        if args.merge:
            header = df_new_rows["ctg_name"][0]
            args.output_fa.write(f">{header}\n")

        for row in df_new_rows.iter_rows(named=True):
            fa_file = query_fa if row["name"] == "qname" else ref_fa
            if not args.merge:
                new_header = f'>{row["ctg_name"]}:{row["st"]}-{row["end"]}\n'
                args.output_fa.write(new_header)
                args.output_fa.write(fa_file.fetch(row["ctg_name"], start=row["st"], end=row["end"]))
                args.output_fa.write("\n")
            else:
                args.output_fa.write(fa_file.fetch(row["ctg_name"], start=row["st"], end=row["end"]))


if __name__ == "__main__":
    raise SystemExit(main())
