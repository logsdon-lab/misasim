use std::{fs::File, io::BufReader};

use clap::Parser;
use eyre::bail;
use itertools::Itertools;
use log::{info, LevelFilter};
use noodles::fasta::{self};
use rand::{rngs::StdRng, seq::SliceRandom, SeedableRng};
use regex::{self, Regex};
use rust_lapper::{Interval, Lapper};
use simple_logger::SimpleLogger;

mod breaks;
mod cli;
mod false_dupe;
mod inversion;
mod io;
mod misassembly;
mod misjoin;
mod sequence;
mod utils;

use crate::{
    io::{write_misassembly_bed, write_new_fasta},
    misassembly::create_all_sequences,
};

use {
    cli::Cli,
    io::{get_outfile_writers, get_regions, Fasta},
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
        .and_then(|f| f.map(BufReader::new).ok());
    let input_regions = get_regions(reader_bed);

    let (output_fa, mut writer_bed) = get_outfile_writers(cli.outfile, cli.outbedfile)?;
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

            let record_interval = Interval {
                start: 1,
                stop: record_length.try_into()?,
                val: (),
            };
            let def_record_regions = Lapper::new(vec![record_interval]);
            let record_regions = input_regions
                .as_ref()
                .and_then(|r| r.get(record_name))
                .unwrap_or(&def_record_regions);

            info!("Processing record: {:?}.", record_name);
            info!("With regions: {:?}.", record_regions);

            let seq = std::str::from_utf8(record.sequence().as_ref())?;

            let seq_segments =
                create_all_sequences(&command, seq, record_regions, seed, randomize_length)?;

            write_new_fasta(record_name, seq, &seq_segments, &mut writer_fa)?;

            if let Some(writer_bed) = writer_bed.as_mut() {
                write_misassembly_bed(record_name, &seq_segments, writer_bed)?;
            }
        }
    }
    Ok(())
}

fn main() -> eyre::Result<()> {
    SimpleLogger::new().with_level(LevelFilter::Debug).init()?;
    let cli = Cli::parse();
    info!("Running the following command:\n{:#?}", cli.command);

    generate_misassemblies(cli)?;
    info!("Completed generating misassemblies.");
    Ok(())
}
