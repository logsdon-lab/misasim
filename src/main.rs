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
mod false_dupe;
mod misjoin;
mod utils;

use {
    breaks::{generate_breaks, write_breaks},
    cli::Cli,
    false_dupe::generate_false_duplication,
    misjoin::generate_deletion,
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

fn read_regions(
    mut reader_bed: Option<bed::Reader<BufReader<File>>>,
) -> Option<HashMap<String, IntervalSet<Position>>> {
    reader_bed.as_mut().map(|input_bed| {
        let mut regions: HashMap<String, IntervalSet<Position>> = HashMap::new();
        for rec in input_bed.records::<3>().flatten() {
            let region = rec.start_position()..rec.end_position();
            regions
                .entry(rec.to_string())
                .and_modify(|r| {
                    r.insert(region);
                })
                .or_default();
        }
        regions
    })
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
    let regions = read_regions(input_bed);

    let mut writer_fa = fasta::Writer::new(output_fa);
    // TODO: async for concurrent record reading.
    for record in reader_fa.records().flatten() {
        let record_name = std::str::from_utf8(record.definition().name())?;
        let record_interval =
            Position::new(1).unwrap()..Position::new(record.sequence().len()).unwrap();
        let def_record_regions = IntervalSet::from_iter(std::iter::once(record_interval));
        let record_regions = regions
            .as_ref()
            .and_then(|r| r.get(record_name))
            .unwrap_or(&def_record_regions);

        info!("Processing record: {:?}.", record_name);
        info!("With regions: {:?}.", record_regions);

        let seq = std::str::from_utf8(record.sequence().as_ref())?;

        match command {
            cli::Commands::Misjoin { number, length, .. }
            | cli::Commands::Gap { number, length, .. } => {
                let deleted_seq = generate_deletion(
                    seq,
                    record_regions,
                    length,
                    number,
                    // If gap, mask deletion.
                    std::mem::discriminant(&command)
                        == std::mem::discriminant(&cli::Commands::Gap {
                            number,
                            length,
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
            cli::Commands::FalseDuplication {
                number,
                length,
                max_duplications,
                seed,
                ..
            } => {
                let false_dupe_seq = generate_false_duplication(
                    seq,
                    record_regions,
                    length,
                    number,
                    max_duplications,
                    seed,
                )?;
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
                let seq_breaks = generate_breaks(seq, record_regions, number, seed)?;
                write_breaks(record_name, seq_breaks, &mut writer_fa, &mut output_bed)?;
            }
        }
    }

    Ok(())
}

fn main() -> eyre::Result<()> {
    SimpleLogger::new().with_level(LevelFilter::Debug).init()?;
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
