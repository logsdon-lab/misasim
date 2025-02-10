use std::{fs::File, io::BufReader};

use clap::Parser;
use eyre::bail;
use iset::IntervalSet;
use itertools::Itertools;
use log::{info, LevelFilter};
use noodles::{
    bed,
    core::Position,
    fasta::{self},
};
use rand::{rngs::StdRng, seq::SliceRandom, SeedableRng};
use regex::{self, Regex};
use simple_logger::SimpleLogger;

mod breaks;
mod cli;
mod false_dupe;
mod inversion;
mod io;
mod misjoin;
mod utils;

use {
    breaks::{generate_breaks, write_breaks},
    cli::Cli,
    false_dupe::generate_false_duplication,
    inversion::generate_inversion,
    io::{get_outfile_writers, get_regions, Fasta},
    misjoin::generate_deletion,
    utils::write_misassembly,
};

fn generate_misassemblies(cli: cli::Cli) -> eyre::Result<()> {
    let command = cli.command;

    let Some(infile) = cli.infile else {
        bail!("No input fasta provided.")
    };
    let mut reader_fa = Fasta::new(infile)?;

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

    let record_groups = reader_fa.lengths();

    let rgx = cli
        .group_by
        .as_deref()
        .map(|rgx| Regex::new(rgx).unwrap())
        .unwrap_or_else(|| Regex::new(".*?").unwrap());

    // Group names by captured groups.
    // ex. [chr10_mat, chr10_pat]
    // * "^.*?_(?<hap>.*?)$" with group by haplotype.
    // * "^(?<chr>.*?)_.*?$" will group by chromosome.
    // * ".*?" will not group as all groups are unique.
    let groups = record_groups
        .into_iter()
        // Sort first by name.
        .sorted_by(|a, b| a.0.cmp(&b.0))
        .chunk_by(|(rec, _)| {
            rgx.captures(rec).map(|captures| {
                captures
                    .iter()
                    .enumerate()
                    .flat_map(|(i, cap)| {
                        // Skip entire string match.
                        if i == 0 {
                            None
                        } else {
                            cap.map(|c| c.as_str().to_owned())
                        }
                    })
                    .collect_vec()
            })
        });

    let mut rng = seed.map_or(StdRng::from_entropy(), StdRng::seed_from_u64);
    for (grp, grps) in &groups {
        if cli.group_by.is_some() {
            log::info!("Grouping by: {grp:?}")
        }
        let grps = grps.collect_vec();
        // Choose one record per group to generate misassemblies.
        let Some(misasm_rec) = grps.choose(&mut rng) else {
            continue;
        };
        for rec in grps.iter() {
            let record_name = &rec.0;
            let record_length: u32 = rec.1.try_into()?;
            let record = reader_fa.fetch(record_name, 1, record_length)?;

            // If not chosen misassembled sequence, then just write record as is.
            if rec != misasm_rec {
                writer_fa.write_record(&record)?;
                continue;
            }

            let record_interval =
                Position::new(1).unwrap()..Position::new(record_length.try_into()?).unwrap();
            let def_record_regions = IntervalSet::from_iter(std::iter::once(record_interval));
            let record_regions = input_regions
                .as_ref()
                .and_then(|r| r.get(record_name))
                .unwrap_or(&def_record_regions);

            info!("Processing record: {:?}.", record_name);
            info!("With regions: {:?}.", record_regions);

            let seq = std::str::from_utf8(record.sequence().as_ref())?;

            match command {
                cli::Commands::Misjoin { number, length }
                | cli::Commands::Gap { number, length } => {
                    let is_gap = std::mem::discriminant(&command)
                        == std::mem::discriminant(&cli::Commands::Gap { number, length });
                    let deleted_seq = generate_deletion(
                        seq,
                        record_regions,
                        length,
                        number,
                        // If gap, mask deletion.
                        is_gap,
                        seed,
                        randomize_length,
                    )?;
                    info!("{} sequence(s) removed.", deleted_seq.removed_seqs.len());

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
                        "{} sequence(s) duplicated.",
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
                cli::Commands::Inversion { number, length } => {
                    let inverted_seq = generate_inversion(
                        seq,
                        record_regions,
                        length,
                        number,
                        seed,
                        randomize_length,
                    )?;
                    info!("{} sequence(s) inverted.", inverted_seq.inverted_seqs.len());

                    write_misassembly(
                        inverted_seq.seq.into_bytes(),
                        inverted_seq.inverted_seqs,
                        record.definition().clone(),
                        &mut writer_fa,
                        output_bed.as_mut(),
                    )?;
                } // cli::Commands::Multiple { subcommands } => {
                  //     unimplemented!()
                  // }
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
