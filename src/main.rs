use clap::{CommandFactory, Parser};
use log::{info, LevelFilter};
use noodles::{
    bed::{self, record::Builder},
    fasta::{self as fasta, reader::Records, record::Sequence},
};
use simple_logger::SimpleLogger;
use std::{
    fs::File,
    io::{stdin, stdout, BufRead, BufReader, IsTerminal, Write},
    path::PathBuf,
};
mod cli;
mod collapse;
mod misjoin;

use {
    cli::{Cli, Commands},
    collapse::{find_all_repeats, generate_collapse},
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

                writer_fa.write_record(&fasta::Record::new(
                    record.definition().clone(),
                    Sequence::from(deleted_seq.seq.into_bytes()),
                ))?;

                let Some(writer_bed) = &mut output_bed else {
                    continue;
                };
                for del_seq in deleted_seq.removed_seqs {
                    // Add deleted sequence to BED file.
                    let record = TryInto::<Builder<3>>::try_into(del_seq)?
                        .set_reference_sequence_name(record_name)
                        .build()?;
                    writer_bed.write_record(&record)?;
                }
            }
            cli::Commands::Collapse {
                length,
                num_repeats,
            } => {
                let repeats = find_all_repeats(seq, length);
                info!("{} repeats found.", repeats.len());

                let collapsed_seq = generate_collapse(seq, &repeats, num_repeats, seed)?;
                info!("{} repeats collapsed.", collapsed_seq.repeats.len());

                writer_fa.write_record(&fasta::Record::new(
                    record.definition().clone(),
                    Sequence::from(collapsed_seq.seq.into_bytes()),
                ))?;

                // Write the BED file if provided.
                let Some(writer_bed) = &mut output_bed else {
                    continue;
                };
                for rp in collapsed_seq.repeats {
                    let record = Into::<Builder<3>>::into(rp)
                        .set_reference_sequence_name(record_name)
                        .build()?;
                    writer_bed.write_record(&record)?;
                }
            }
            cli::Commands::FalseDuplication { length, number } => {
                use rand::prelude::*;
                let mut rng = seed.map_or(StdRng::from_entropy(), StdRng::seed_from_u64);
                let mut false_dupe_pos = (0..seq.len()).choose_multiple(&mut rng, number);
                false_dupe_pos.sort();

                println!("False duplication with length: {}", length);
            }
            cli::Commands::Break { number } => {}
        }
    }
    Ok(())
}

fn main() -> eyre::Result<()> {
    SimpleLogger::new().with_level(LevelFilter::Info).init()?;

    let (file, out_fa, out_bed, cmd, seed) =
        if std::env::var("DEBUG").map_or(false, |v| v == "1" || v == "true") {
            let cmd = Commands::Misjoin { number: 10 };
            let file = PathBuf::from("test/data/HG00171_chr9_haplotype2-0000142.fa");
            let out_fa: Option<PathBuf> = Some(PathBuf::from("output.fa"));
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
    Ok(())
}
