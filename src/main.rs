use bio::io::fasta::{FastaRead, Reader, Records};
use clap::{CommandFactory, Parser};
use std::io::{stdin, BufRead, BufReader};
use std::path::PathBuf;
use std::{fs::File, io::IsTerminal};
use suffix::SuffixTable;

mod cli;

use cli::{Cli, Commands};

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
        let curr_sfx_len = curr_sfx.len();

        if let Some((sfx, sfx_pos)) = prev_sfx {

            // Does current suffix start with this suffix?
            if !curr_sfx.starts_with(sfx) {
                prev_sfx = Some((curr_sfx, 0));

                if let Some(repeat) = prev_repeat {
                    repeats.push(repeat);
                }
                prev_repeat = None;
                continue;
            }
            let prev_sfx_len = sfx.len();
            // Check slice ahead.
            let sfx_end = sfx_pos + prev_sfx_len;
            let Some(next_sfx) = &curr_sfx.get(sfx_end..sfx_end + prev_sfx_len) else {
                prev_sfx = None;
                continue;
            };

            if *next_sfx == sfx {
                if let Some(repeat) = prev_repeat.as_mut() {
                    repeat.count += 1;
                } else {
                    prev_repeat = Some(Repeat {
                        seq: sfx,
                        // Calculate the start of the repeat. If overflows, means idx is already at the start.
                        start: (*idx as usize).checked_sub(prev_sfx_len).unwrap_or(*idx as usize),
                        count: 2,
                    });
                }
            } else {
                continue;
            }
            prev_sfx = Some((sfx, sfx_end));

        } else {
            if curr_sfx_len < min_len {
                continue;
            }
            prev_sfx = Some((curr_sfx, 0));
        }
    }
    repeats
}

fn generate_misjoin(seq: &str) -> eyre::Result<String> {
    let mut new_seq = String::new();

    Ok(new_seq)
}

fn generate_collapse(
    seq: &str,
    repeats: &[Repeat],
    num_repeats: usize,
) -> eyre::Result<String> {
    let mut new_seq = String::new();

    for repeat in repeats {}
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
        let seq = std::str::from_utf8(record.seq())?;
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
            },
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
            num_repeats: 2,
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
        let fasta_reader = Reader::new(buf_reader);
        let results = generate_misassemblies(fasta_reader.records(), cmd);
    } else {
        let buf_reader = BufReader::new(File::open(&file).unwrap());
        let fasta_reader = Reader::new(buf_reader);
        let results = generate_misassemblies(fasta_reader.records(), cmd);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_repeats() {
        let seq = "ATTTTATTTTAATTTTAATTTTAATTTT";
        let stbl = SuffixTable::new(seq);
        println!("{stbl:?}");

        let repeats = find_repeats(&stbl, 5);
        assert!(repeats.len() == 2);
    }
}
