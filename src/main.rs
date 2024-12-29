use std::{fs::File, io::BufReader};

use clap::Parser;
use eyre::bail;
use iset::IntervalSet;
use log::{info, LevelFilter};
use noodles::{bed, core::Position, fasta};
use simple_logger::SimpleLogger;

mod breaks;
mod cli;
mod false_dupe;
mod io;
mod misjoin;
mod utils;

use {
    breaks::{generate_breaks, write_breaks},
    cli::Cli,
    false_dupe::generate_false_duplication,
    io::{get_fa_reader, get_outfile_writers, get_regions},
    misjoin::generate_deletion,
    utils::write_misassembly,
};

fn generate_misassemblies(cli: cli::Cli) -> eyre::Result<()> {
    let command = cli.command;

    let Some(infile) = cli.infile else {
        bail!("No input fasta provided.")
    };
    let mut reader_fa = fasta::Reader::new(get_fa_reader(infile)?);

    // https://rust-cli.github.io/book/in-depth/machine-communication.html
    let reader_bed = cli
        .inbedfile
        .as_ref()
        .map(File::open)
        .and_then(|f| f.map(BufReader::new).ok())
        .map(bed::Reader::new);
    let input_regions = get_regions(reader_bed);

    let (output_fa, mut output_bed) = get_outfile_writers(cli.outfile, cli.outbedfile)?;
    let mut writer_fa = fasta::Writer::new(output_fa);

    let seed = cli.seed;
    let randomize_length = cli.randomize_length;
    if let Some(seed) = seed {
        log::info!("Random seed: {seed:?}");
    } else {
        log::info!("No random seed provided. Generating a random seed per record.");
    }
    log::info!("Randomizing length: {randomize_length}");

    // TODO: async for concurrent record reading.
    for record in reader_fa.records().flatten() {
        let record_name = std::str::from_utf8(record.definition().name())?;
        let record_interval =
            Position::new(1).unwrap()..Position::new(record.sequence().len()).unwrap();
        let def_record_regions = IntervalSet::from_iter(std::iter::once(record_interval));
        let record_regions = input_regions
            .as_ref()
            .and_then(|r| r.get(record_name))
            .unwrap_or(&def_record_regions);

        info!("Processing record: {:?}.", record_name);
        info!("With regions: {:?}.", record_regions);

        let seq = std::str::from_utf8(record.sequence().as_ref())?;

        match command {
            cli::Commands::Misjoin { number, length } | cli::Commands::Gap { number, length } => {
                let deleted_seq = generate_deletion(
                    seq,
                    record_regions,
                    length,
                    number,
                    // If gap, mask deletion.
                    std::mem::discriminant(&command)
                        == std::mem::discriminant(&cli::Commands::Gap { number, length }),
                    seed,
                    randomize_length,
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
            } => {
                let false_dupe_seq = generate_false_duplication(
                    seq,
                    record_regions,
                    length,
                    number,
                    max_duplications,
                    seed,
                    randomize_length,
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
            cli::Commands::Break { number, .. } => {
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
    // let cli = if std::env::var("DEBUG").map_or(false, |v| v == "1" || v == "true") {
    //     Cli {
    //         command: Commands::Break {
    //             number: 10,
    //             length: 5000,
    //             infile: PathBuf::from("test/data/HG00171_chr9_haplotype2-0000142_full.fa"),
    //             inbedfile: Some(PathBuf::from("test/data/region.bed")),
    //             outfile: None,
    //             outbedfile: None,
    //             seed: Some(42),
    //         },
    //     }
    // } else {
    //     Cli::parse()
    // };
    info!("Running the following command:\n{:#?}", cli.command);

    generate_misassemblies(cli)?;
    info!("Completed generating misassemblies.");
    Ok(())
}
