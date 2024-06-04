use clap::{CommandFactory, Parser};
use itertools::Itertools;
use log::{info, LevelFilter};
use noodles::{
    bed::{
        self,
        record::{Builder, OptionalFields},
    },
    core::Position,
    fasta::{self as fasta, reader::Records, record::Definition},
};
use simple_logger::SimpleLogger;
use std::{
    fs::File,
    io::{stdin, stdout, BufRead, BufReader, IsTerminal, Write},
    path::PathBuf,
};
use utils::{get_sequence_segments, write_misassembly};
mod cli;
mod collapse;
mod false_dupe;
mod misjoin;
mod utils;

use {
    cli::{Cli, Commands},
    collapse::generate_collapse,
    false_dupe::generate_false_duplication,
    misjoin::generate_deletion,
};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

fn generate_misassemblies<B: BufRead, O: Write>(
    records: Records<B>,
    output_fa: O,
    mut output_bed: Option<bed::Writer<File>>,
    seed: Option<u64>,
    command: cli::Commands,
) -> eyre::Result<()> {
    let mut writer_fa = fasta::Writer::new(output_fa);

    for record in records.into_iter() {
        let record = record?;
        let full_record_name = record.definition().to_string();
        let record_name = full_record_name
            .strip_prefix('>')
            .unwrap_or(&full_record_name);

        info!("Processing record: {}.", record_name);

        let seq = std::str::from_utf8(record.sequence().as_ref())?;
        match command {
            cli::Commands::Misjoin { number } | cli::Commands::Gap { number } => {
                let deleted_seq = generate_deletion(
                    seq,
                    number,
                    // If gap, mask deletion.
                    command == cli::Commands::Gap { number },
                    seed,
                )?;
                write_misassembly(
                    deleted_seq.seq.into_bytes(),
                    deleted_seq
                        .removed_seqs
                        .into_iter()
                        .flat_map(TryInto::<Builder<3>>::try_into),
                    record.definition().clone(),
                    &mut writer_fa,
                    output_bed.as_mut(),
                )?;
            }
            cli::Commands::Collapse {
                length,
                num_repeats,
            } => {
                let collapsed_seq = generate_collapse(seq, length, num_repeats, seed)?;
                info!("{} repeats collapsed.", collapsed_seq.repeats.len());

                write_misassembly(
                    collapsed_seq.seq.into_bytes(),
                    collapsed_seq
                        .repeats
                        .into_iter()
                        .flat_map(TryInto::<Builder<3>>::try_into),
                    record.definition().clone(),
                    &mut writer_fa,
                    output_bed.as_mut(),
                )?;
            }
            cli::Commands::FalseDuplication {
                number,
                max_duplications,
            } => {
                let false_dupe_seq =
                    generate_false_duplication(seq, number, max_duplications, seed)?;

                write_misassembly(
                    false_dupe_seq.seq.into_bytes(),
                    false_dupe_seq
                        .duplicated_seqs
                        .into_iter()
                        .flat_map(TryInto::<Builder<3>>::try_into),
                    record.definition().clone(),
                    &mut writer_fa,
                    output_bed.as_mut(),
                )?;
            }
            cli::Commands::Break { number } => {
                let seq_segments = get_sequence_segments(seq, number, seed)
                    .unwrap()
                    .collect_vec();
                for (i, (start, end, _)) in seq_segments.into_iter().enumerate() {
                    let bseq = seq[start..end].bytes().collect_vec();
                    let region = bed::Record::<3>::builder()
                        .set_start_position(Position::new(start).unwrap())
                        .set_optional_fields(OptionalFields::from(vec!["Break".to_string()]));
                    let new_definiton = format!("{}_{}", record.definition(), i);
                    write_misassembly(
                        bseq,
                        std::iter::once(region),
                        Definition::new(new_definiton, None),
                        &mut writer_fa,
                        output_bed.as_mut(),
                    )?;
                }
            }
        }
    }
    Ok(())
}

fn main() -> eyre::Result<()> {
    SimpleLogger::new().with_level(LevelFilter::Info).init()?;

    let (file, out_fa, out_bed, cmd, seed) =
        if std::env::var("DEBUG").map_or(false, |v| v == "1" || v == "true") {
            let cmd = Commands::Break { number: 10 };
            let file = PathBuf::from("test/data/HG00171_chr9_haplotype2-0000142.fa");
            let out_fa: Option<PathBuf> = Some(PathBuf::from("test/output/output.fa"));
            let out_bed = Some(PathBuf::from("test/data/test.bed"));
            let seed = Some(42);
            (file, out_fa, out_bed, cmd, seed)
        } else {
            let cli = Cli::parse();
            let file = cli.infile;
            let out_fa = cli.outfile;
            let out_bed = cli.outbedfile;
            let cmd = cli.command;
            let seed = cli.seed;
            (file, out_fa, out_bed, cmd, seed)
        };

    info!("Using the following configuration:");
    info!("Command: {cmd:?}.");
    info!("Seed: {seed:?}.");

    // TODO: This won't work with rayon.
    let output_fa: Box<dyn Write> = if let Some(outfile) = out_fa {
        Box::new(File::create(outfile)?)
    } else {
        Box::new(stdout().lock())
    };
    let output_bed = out_bed
        .and_then(|f| File::create(f).ok())
        .map(bed::Writer::new);

    // https://rust-cli.github.io/book/in-depth/machine-communication.html
    if file == PathBuf::from("-") {
        if stdin().is_terminal() {
            Cli::command().print_help()?;
            std::process::exit(2);
        }
        info!("Reading from stdin.");

        let buf_reader = BufReader::new(stdin().lock());
        let mut fasta_reader = fasta::Reader::new(buf_reader);
        generate_misassemblies(fasta_reader.records(), output_fa, output_bed, seed, cmd)?;
    } else {
        info!("Reading {file:?}.");

        let buf_reader = BufReader::new(File::open(&file).unwrap());
        let mut fasta_reader = fasta::Reader::new(buf_reader);
        generate_misassemblies(fasta_reader.records(), output_fa, output_bed, seed, cmd)?;
    }
    info!("Completed generating misassemblies.");
    Ok(())
}
