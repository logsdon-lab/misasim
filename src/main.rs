use breaks::write_breaks;
use clap::{CommandFactory, Parser};
use iset::IntervalSet;
use log::{info, LevelFilter};
use noodles::{bed, core::Position, fasta};
use simple_logger::SimpleLogger;
use std::{
    collections::HashMap,
    fs::File,
    io::{stdin, stdout, BufRead, BufReader, IsTerminal, Write},
    path::PathBuf,
};

use utils::write_misassembly;
mod breaks;
mod cli;
mod collapse;
mod false_dupe;
mod misjoin;
mod utils;

use {
    breaks::generate_breaks, cli::Cli, collapse::generate_collapse,
    false_dupe::generate_false_duplication, misjoin::generate_deletion,
};

type Outfiles = (Box<dyn Write>, Option<bed::Writer<File>>);

fn get_outfiles(outfile: Option<PathBuf>, outbedfile: Option<PathBuf>) -> eyre::Result<Outfiles> {
    // TODO: This won't work with rayon.
    let output_fa: Box<dyn Write> = if let Some(outfile) = outfile {
        Box::new(File::create(outfile)?)
    } else {
        Box::new(stdout().lock())
    };
    let output_bed = outbedfile
        .and_then(|f| File::create(f).ok())
        .map(bed::Writer::new);

    Ok((output_fa, output_bed))
}

fn get_records(infile: PathBuf) -> eyre::Result<Box<dyn BufRead>> {
    if infile == PathBuf::from("-") {
        if stdin().is_terminal() {
            Cli::command().print_help()?;
            std::process::exit(2);
        }
        info!("Reading from stdin.");
        Ok(Box::new(BufReader::new(stdin().lock())) as Box<dyn BufRead>)
    } else {
        info!("Reading {infile:?}.");
        Ok(Box::new(BufReader::new(File::open(&infile).unwrap())))
    }
}

fn generate_misassemblies(command: cli::Commands) -> eyre::Result<()> {
    let (output_fa, mut output_bed, input_records_buf, input_bed, seed) = match &command {
        cli::Commands::Misjoin {
            infile,
            inbedfile,
            outfile,
            outbedfile,
            seed,
            ..
        }
        | cli::Commands::Gap {
            infile,
            inbedfile,
            outfile,
            outbedfile,
            seed,
            ..
        }
        | cli::Commands::Collapse {
            infile,
            inbedfile,
            outfile,
            outbedfile,
            seed,
            ..
        }
        | cli::Commands::FalseDuplication {
            infile,
            inbedfile,
            outfile,
            outbedfile,
            seed,
            ..
        }
        | cli::Commands::Break {
            infile,
            inbedfile,
            outfile,
            outbedfile,
            seed,
            ..
        } => {
            let (output_fa, output_bed) = get_outfiles(outfile.clone(), outbedfile.clone())?;
            // https://rust-cli.github.io/book/in-depth/machine-communication.html
            let buf = get_records(infile.clone())?;
            let input_bed = inbedfile
                .as_ref()
                .map(File::open)
                .and_then(|f| f.map(BufReader::new).ok())
                .map(bed::Reader::new);
            (output_fa, output_bed, buf, input_bed, seed)
        }
    };
    let mut reader_fa = fasta::Reader::new(input_records_buf);
    if let Some(mut reader_bed) = input_bed {
        let mut regions: HashMap<String, IntervalSet<Position>> = HashMap::new();
        for rec in reader_bed.records::<3>().flatten() {
            let region = rec.start_position()..rec.end_position();
            regions
                .entry(rec.to_string())
                .and_modify(|r| {
                    r.insert(region);
                })
                .or_default();
        }
    }

    let mut writer_fa = fasta::Writer::new(output_fa);
    // TODO: async for concurrent record reading.
    for record in reader_fa.records().flatten() {
        let record_name = std::str::from_utf8(record.definition().name())?;
        info!("Processing record: {:?}.", record_name);

        let seq = std::str::from_utf8(record.sequence().as_ref())?;

        match command {
            cli::Commands::Misjoin { number, .. } | cli::Commands::Gap { number, .. } => {
                let deleted_seq = generate_deletion(
                    seq,
                    number,
                    // If gap, mask deletion.
                    std::mem::discriminant(&command)
                        == std::mem::discriminant(&cli::Commands::Gap {
                            number,
                            infile: PathBuf::default(),
                            inbedfile: None,
                            outfile: None,
                            outbedfile: None,
                            seed: None,
                        }),
                    *seed,
                )?;
                info!("{} sequences removed.", deleted_seq.removed_seqs.len());

                write_misassembly(
                    deleted_seq.seq.into_bytes(),
                    deleted_seq.removed_seqs,
                    record.definition().clone(),
                    &mut writer_fa,
                    output_bed.as_mut(),
                )?;
            }
            cli::Commands::Collapse {
                length,
                number,
                seed,
                ..
            } => {
                let collapsed_seq = generate_collapse(seq, length, number, seed)?;
                info!("{} repeats collapsed.", collapsed_seq.repeats.len());

                write_misassembly(
                    collapsed_seq.seq.into_bytes(),
                    collapsed_seq.repeats,
                    record.definition().clone(),
                    &mut writer_fa,
                    output_bed.as_mut(),
                )?;
            }
            cli::Commands::FalseDuplication {
                number,
                max_duplications,
                seed,
                ..
            } => {
                let false_dupe_seq =
                    generate_false_duplication(seq, number, max_duplications, seed)?;
                info!(
                    "{} sequences duplicated.",
                    false_dupe_seq.duplicated_seqs.len()
                );

                write_misassembly(
                    false_dupe_seq.seq.into_bytes(),
                    false_dupe_seq.duplicated_seqs,
                    record.definition().clone(),
                    &mut writer_fa,
                    output_bed.as_mut(),
                )?;
            }
            cli::Commands::Break { number, seed, .. } => {
                let seq_breaks = generate_breaks(seq, number, seed)?;
                write_breaks(record_name, seq_breaks, &mut writer_fa, &mut output_bed)?;
            }
        }
    }

    Ok(())
}

fn main() -> eyre::Result<()> {
    SimpleLogger::new().with_level(LevelFilter::Info).init()?;
    let cli = Cli::parse();
    info!("Running the following command:\n{:#?}", cli.command);

    // if std::env::var("DEBUG").map_or(false, |v| v == "1" || v == "true") {
    //     let cmd = Commands::Break { number: 10 };
    //     let file = PathBuf::from("test/data/HG00171_chr9_haplotype2-0000142.fa");
    //     let out_fa: Option<PathBuf> = Some(PathBuf::from("test/output/output.fa"));
    //     let out_bed = Some(PathBuf::from("test/data/test.bed"));
    //     let seed = Some(42);
    // }

    generate_misassemblies(cli.command)?;
    info!("Completed generating misassemblies.");
    Ok(())
}
