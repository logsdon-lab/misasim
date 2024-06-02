use clap::{CommandFactory, Parser};
use iset::IntervalMap;
use itertools::Itertools;
use noodles::fasta::{reader::Records, Reader};

use std::{
    collections::HashSet,
    fs::File,
    io::{stdin, BufRead, BufReader, IsTerminal},
    path::PathBuf,
};

mod cli;

use cli::{Cli, Commands};

#[derive(Debug, PartialEq, Eq, Clone, Hash)]
struct Repeat {
    seq: String,
    start: usize,
    count: usize,
}

fn find_all_repeats(seq: &str, min_len: usize) -> HashSet<Repeat> {
    let rev_seq = seq.chars().rev().collect::<String>();
    let fwd_repeats = find_repeats(seq, min_len);
    let rev_repeats = find_repeats(&rev_seq, min_len);

    // Reorient reverse repeats.
    rev_repeats
        .into_iter()
        .map(|mut repeat| {
            repeat.start = seq.len() - repeat.start - (repeat.seq.len() * repeat.count);
            repeat.seq = repeat.seq.chars().rev().collect();
            repeat
        })
        .chain(fwd_repeats)
        .collect()
}

fn find_repeats(seq: &str, min_len: usize) -> Vec<Repeat> {
    let mut repeats: Vec<Repeat> = vec![];
    let mut prev_sfx: Option<(&str, usize)> = None;
    let mut prev_repeat: Option<Repeat> = None;

    for (idx, curr_sfx) in seq
        .char_indices()
        .flat_map(|(i, _)| seq.get(i..).map(|s| (i, s)))
        .sorted_by(|a, b| a.1.cmp(b.1))
    {
        if let Some((mut sfx, sfx_pos)) = prev_sfx {
            // Does current suffix start with this suffix?
            if !curr_sfx.starts_with(sfx) {
                // Remove prev suffix from current suffix.
                // This isolates repeated elements at start of current suffix.
                prev_sfx = curr_sfx
                    .strip_suffix(sfx)
                    .map_or(Some((curr_sfx, 0)), |s| Some((s, 0)));

                if let Some((s, _)) = prev_sfx {
                    sfx = s;
                }
                if let Some(repeat) = prev_repeat {
                    repeats.push(repeat);
                    prev_repeat = None;
                }
            }
            let prev_sfx_len = sfx.len();
            // Remove small suffixes.
            if prev_sfx_len < min_len {
                prev_sfx = None;
                continue;
            }
            // Check slice ahead.
            let sfx_end = sfx_pos + prev_sfx_len;
            let Some(next_sfx) = &curr_sfx.get(sfx_end..sfx_end + prev_sfx_len) else {
                continue;
            };

            // If same as previous suffix, increment count of repeat.
            // Otherwise, store repeat.
            if *next_sfx == sfx {
                if let Some(repeat) = prev_repeat.as_mut() {
                    repeat.count += 1;
                } else {
                    prev_repeat = Some(Repeat {
                        seq: sfx.to_owned(),
                        start: idx,
                        count: 2,
                    });
                }
            } else {
                continue;
            }
            prev_sfx = Some((sfx, sfx_end));
        } else {
            prev_sfx = Some((curr_sfx, 0));
        }
    }

    if let Some(repeat) = prev_repeat {
        repeats.push(repeat);
    }

    repeats
}

fn generate_misjoin(seq: &str) -> eyre::Result<String> {
    let mut new_seq = String::new();

    Ok(new_seq)
}

fn generate_collapse(
    seq: &str,
    repeats: &HashSet<Repeat>,
    num_repeats: usize,
) -> eyre::Result<String> {
    let mut new_seq: String = String::with_capacity(seq.len());

    let intervals: IntervalMap<usize, Repeat> =
        IntervalMap::from_iter(repeats.iter().map(|repeat| {
            let (start, stop) = (
                repeat.start,
                repeat.start + (repeat.seq.len() * repeat.count),
            );
            (start..stop, repeat.clone())
        }));

    let mut curr_sinterval: Option<core::ops::Range<usize>> = None;
    let mut curr_einterval: Option<core::ops::Range<usize>> = None;
    // Iterate thru intervals to construct slices.
    for (range, repeat) in intervals.iter(0..seq.len()) {
        println!("Current: {repeat:?} {range:?}");
        let (start, end) = (range.start, range.end);

        // Add starting sequence before the repeat.
        for si in intervals
            .overlap(start.saturating_sub(1))
            .filter(|(r, _)| r.start < start && r.end <= start)
        {
            // new_seq.push_str(seq.get(i..start).unwrap_or_default());
            println!("{:?}", si);
        }

        // Collapse the repeat.
        new_seq.push_str(&repeat.seq);

        // Add the remaining repeats if num_repeats less than repeat count.
        let remaining_repeats = seq
            .get(end - (repeat.seq.len() * repeat.count.saturating_sub(num_repeats))..end)
            .unwrap_or_default();

        new_seq.push_str(remaining_repeats);

        // for ei in intervals.overlap(end) {

        // }
        // println!("Interval: {start} - {stop}");
        // if let Some(next_repeat) = sorted_repeats.peek() {
        //     println!("Next: {next_repeat:?}");
        // } else {
        // }
    }
    // println!("{:?}", intervals);

    Ok(new_seq)
}

fn generate_false_duplication(seq: &str) -> eyre::Result<String> {
    let mut new_seq = String::new();

    Ok(new_seq)
}

fn generate_misassemblies<B: BufRead>(
    records: Records<B>,
    command: cli::Commands,
) -> eyre::Result<()> {
    for record in records {
        let record = record?;
        let seq = std::str::from_utf8(record.sequence().as_ref())?;
        match command {
            cli::Commands::Misjoin { length, number } => {
                println!("Misjoin with length: {}", length);
            }
            cli::Commands::Collapse {
                length,
                num_repeats,
            } => {
                let repeats = find_all_repeats(seq, length);
                let seqs = generate_collapse(seq, &repeats, num_repeats);
            }
            cli::Commands::FalseDuplication { length, number } => {
                println!("False duplication with length: {}", length);
            }
            cli::Commands::Gap { length, number } => {
                println!("Gap with length: {}", length);
            }
        }
    }

    Ok(())
}

fn main() -> eyre::Result<()> {
    let (file, cmd) = if std::env::var("DEBUG").map_or(false, |v| v == "1" || v == "true") {
        let cmd = Commands::Collapse {
            length: 5,
            num_repeats: 1,
        };
        let file = PathBuf::from("test/data/test.fa");
        (file, cmd)
    } else {
        let cli = Cli::parse();
        let file = cli.infile;
        let cmd = cli.command;
        (file, cmd)
    };

    // https://rust-cli.github.io/book/in-depth/machine-communication.html
    if file == PathBuf::from("-") {
        if stdin().is_terminal() {
            Cli::command().print_help()?;
            std::process::exit(2);
        }

        let buf_reader = BufReader::new(stdin().lock());
        let mut fasta_reader = Reader::new(buf_reader);
        let results = generate_misassemblies(fasta_reader.records(), cmd);
    } else {
        let buf_reader = BufReader::new(File::open(&file).unwrap());
        let mut fasta_reader = Reader::new(buf_reader);
        let results = generate_misassemblies(fasta_reader.records(), cmd);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_repeats() {
        let seq = "ATTTTATTTT";
        let repeats = find_all_repeats(&seq, 5);
        assert_eq!(
            vec![Repeat {
                seq: String::from("ATTTT"),
                start: 0,
                count: 2
            }],
            repeats.into_iter().collect_vec()
        );
    }

    #[test]
    fn test_find_repeats_overlap() {
        let seq = "ATTTTATTTTA";
        let repeats = find_all_repeats(&seq, 5);
        assert_eq!(
            vec![
                Repeat {
                    seq: "TTTTA".to_string(),
                    start: 1,
                    count: 2
                },
                Repeat {
                    seq: "ATTTT".to_string(),
                    start: 0,
                    count: 2
                }
            ],
            repeats.into_iter().collect_vec()
        );
    }

    #[test]
    fn test_find_repeats_multiple() {
        let seq = "GCCCCGCCCCAATTTTAATTTTAATTTT";
        let repeats = find_all_repeats(&seq, 5);
        assert_eq!(
            vec![
                Repeat {
                    seq: "TAATTT".to_string(),
                    start: 15,
                    count: 2
                },
                Repeat {
                    seq: "AATTTT".to_string(),
                    start: 4,
                    count: 3
                },
                Repeat {
                    seq: "AATTTT".to_string(),
                    start: 16,
                    count: 3
                },
                Repeat {
                    seq: "GCCCC".to_string(),
                    start: 0,
                    count: 2
                },
                Repeat {
                    seq: "TTAATT".to_string(),
                    start: 14,
                    count: 2
                },
                Repeat {
                    seq: "TTTTAA".to_string(),
                    start: 12,
                    count: 2
                },
                Repeat {
                    seq: "ATTTTA".to_string(),
                    start: 11,
                    count: 2
                },
                Repeat {
                    seq: "TTTAAT".to_string(),
                    start: 13,
                    count: 2
                }
            ],
            repeats.into_iter().collect_vec()
        )
    }

    #[test]
    fn test_generate_collapse() {
        let seq = "ATTTTATTTT";
        let repeats = find_all_repeats(&seq, 5);
        let new_seq = generate_collapse(seq, &repeats, 20).unwrap();

        assert_eq!("ATTTT", new_seq);
    }

    #[test]
    fn test_generate_collapse_multiple() {
        let seq = "AAAGGCCCGGCCCGGGGATTTTATTTT";
        let repeats = find_all_repeats(&seq, 5);
        let new_seq = generate_collapse(seq, &repeats, 20).unwrap();

        unimplemented!()
    }

    #[test]
    fn test_generate_collapse_interm_seq() {
        let seq = "ATTTTATTTTGCCGAATTTTAATTTTAATTTT";
        let repeats = find_all_repeats(&seq, 5);
        let new_seq = generate_collapse(seq, &repeats, 20).unwrap();
        unimplemented!()
    }
}
