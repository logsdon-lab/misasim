use itertools::Itertools;
use suffix::SuffixTable;
use clap::{CommandFactory, Parser};
use noodles::fasta::{reader::Records, Reader};

use std::{
    path::PathBuf,
    fs::File, io::IsTerminal,
    io::{stdin, BufRead, BufReader}
};

mod cli;

use cli::{Cli, Commands};

#[derive(Debug, PartialEq, Eq)]
struct Repeat<'a> {
    seq: &'a str,
    start: usize,
    count: usize,
}

fn find_repeats<'a>(stbl: &'a SuffixTable, min_len: usize) -> Vec<Repeat<'a>> {
    let mut repeats: Vec<Repeat> = vec![];
    let mut prev_sfx: Option<(&str, usize)> = None;
    let mut prev_repeat: Option<Repeat> = None;

    for (num, idx) in stbl.table().iter().enumerate() {
        let curr_sfx = stbl.suffix(num);

        if let Some((sfx, sfx_pos)) = prev_sfx {
            // Does current suffix start with this suffix?
            if !curr_sfx.starts_with(sfx) {
                // Remove prev suffix from current suffix.
                // This isolates repeated elements
                prev_sfx = curr_sfx.strip_suffix(sfx).map_or(Some((curr_sfx, 0)), |s| Some((s, 0)));

                if let Some(repeat) = prev_repeat {
                    repeats.push(repeat);
                }
                prev_repeat = None;
                continue;
            }
            let prev_sfx_len = sfx.len();
            // Remove small suffixes.
            if prev_sfx_len < min_len {
                prev_sfx = curr_sfx.strip_suffix(sfx).map_or(Some((curr_sfx, 0)), |s| Some((s, 0)));
                continue;
            }
            // Check slice ahead.
            let sfx_end = sfx_pos + prev_sfx_len;
            let Some(next_sfx) = &curr_sfx.get(sfx_end..sfx_end + prev_sfx_len) else {
                prev_sfx = None;
                continue;
            };

            if *next_sfx == sfx {
                // TODO: Incorrect count with intermediate sequence.
                if let Some(repeat) = prev_repeat.as_mut() {
                    repeat.count += 1;
                } else {
                    prev_repeat = Some(Repeat {
                        seq: sfx,
                        // Calculate the start of the repeat. If overflows, means idx is already at the start.
                        start: (*idx as usize)
                            .checked_sub(prev_sfx_len)
                            .unwrap_or(*idx as usize),
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

fn generate_collapse(seq: &str, repeats: &[Repeat], num_repeats: usize) -> eyre::Result<String> {
    let mut new_seq: String = String::with_capacity(seq.len());

    let mut sorted_repeats = repeats
        .iter()
        .sorted_by(|a, b| a.start.cmp(&b.start))
        .peekable();

    // Iterate thru intervals to construct slices.
    let mut i = 0;
    while let Some(repeat) = sorted_repeats.next() {
        println!("Current: {repeat:?}");
        let (start, stop) = (
            repeat.start,
            repeat.start + (repeat.seq.len() * repeat.count),
        );
        // Add starting sequence before the repeat.
        new_seq.push_str(seq.get(i..start).unwrap_or_default());

        // Collapse the repeat.
        new_seq.push_str(repeat.seq);

        // Add the remaining repeats if num_repeats less than repeat count.
        let remaining_repeats = seq
            .get(stop - (repeat.seq.len() * repeat.count.saturating_sub(num_repeats))..stop)
            .unwrap_or_default();

        new_seq.push_str(remaining_repeats);

        println!("Interval: {start} - {stop}");
        if let Some(next_repeat) = sorted_repeats.peek() {
            println!("Next: {next_repeat:?}");
        } else {
        }

        i = stop;
    }

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
                let stbl = SuffixTable::new(seq);
                let repeats = find_repeats(&stbl, length);
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
        let stbl = SuffixTable::new(seq);
        let repeats = find_repeats(&stbl, 5);
        assert_eq!(
            [Repeat {
                seq: "ATTTT",
                start: 0,
                count: 2
            }],
            repeats.as_slice()
        );
    }

    #[test]
    fn test_find_repeats_overlap() {
        let seq = "ATTTTATTTTA";
        let stbl = SuffixTable::new(seq);
        let repeats = find_repeats(&stbl, 5);
        assert_eq!(
            [
                Repeat {
                    seq: "ATTTT",
                    start: 0,
                    count: 2
                },
                Repeat {
                    seq: "TTTTA",
                    start: 1,
                    count: 2
                }
            ],
            repeats.as_slice(),
        )
    }

    #[test]
    fn test_find_repeats_multiple() {
        let seq = "ATTTTATTTTAATTTTAATTTTAATTTT";
        let stbl = SuffixTable::new(seq);
        let repeats = find_repeats(&stbl, 5);
        assert_eq!(
            [
                Repeat {
                    seq: "AATTTT",
                    start: 10,
                    count: 3
                },
                Repeat {
                    seq: "ATTTT",
                    start: 0,
                    count: 2
                },
                Repeat {
                    seq: "TAATTT",
                    start: 9,
                    count: 3
                },
                Repeat {
                    seq: "TTAATT",
                    start: 8,
                    count: 3
                },
                Repeat {
                    seq: "TTTAAT",
                    start: 7,
                    count: 3
                },
                Repeat {
                    seq: "TTTTAA",
                    start: 6,
                    count: 3
                }
            ],
            repeats.as_slice()
        );
    }

    #[test]
    fn test_generate_collapse() {
        let seq = "ATTTTATTTT";
        let stbl = SuffixTable::new(seq);
        let repeats = find_repeats(&stbl, 5);
        let new_seq = generate_collapse(seq, &repeats, 20).unwrap();

        assert_eq!("ATTTT", new_seq);
    }

    #[test]
    fn test_generate_collapse_multiple() {
        let seq = "ATTTTATTTTAATTTTAATTTTAATTTT";
        let stbl = SuffixTable::new(seq);
        let repeats = find_repeats(&stbl, 5);
        let new_seq = generate_collapse(seq, &repeats, 20).unwrap();

        unimplemented!()
    }

    #[test]
    fn test_generate_collapse_interm_seq() {
        let seq = "ATTTTATTTTGCCGAATTTTAATTTTAATTTT";
        let stbl = SuffixTable::new(seq);
        let repeats = find_repeats(&stbl, 5);
        let new_seq = generate_collapse(seq, &repeats, 20).unwrap();
        unimplemented!()
    }
}
