use std::{fs::File, io::IsTerminal};
use std::path::PathBuf;
use std::io::{stdin, BufRead, BufReader};
use clap::{CommandFactory, Parser};
use bio::io::{fastq, fasta};
use bio::io::fastq::FastqRead;
use bio::io::fasta::FastaRead;
use eyre::{bail, Error};
use suffix::SuffixTable;

mod cli;

use cli::Cli;

fn main() -> eyre::Result<()> {
    let cli = Cli::parse();

    match &cli.command {
        Some(command) => match command {
            cli::Commands::Misjoin { length } => {
                println!("Misjoin with length: {}", length);
            },
            cli::Commands::Collapse => todo!(),
            cli::Commands::FalseDuplication { length } => {
                println!("False duplication with length: {}", length);
            },
        },
        None => {
            bail!("No command provided.");
        }
    }

    let file = cli.infile;


    // https://rust-cli.github.io/book/in-depth/machine-communication.html
    if file == PathBuf::from("-") {
        if stdin().is_terminal() {
            Cli::command().print_help()?;
            std::process::exit(2);
        }

        let buf_reader = BufReader::new(stdin().lock());
        let fasta_reader = fasta::Reader::new(buf_reader);
        for record in fasta_reader.records() {
            let record = record?;
            let seq = std::str::from_utf8(record.seq())?;
            println!("ID: {}", record.id());
            println!("Seq: {}", seq);

            let stbl = SuffixTable::new(seq);
        }
    } else {
        let buf_reader = BufReader::new(File::open(&file).unwrap());
        let fasta_reader = fasta::Reader::new(buf_reader);
        for record in fasta_reader.records() {
            let record = record?;
            let seq = std::str::from_utf8(record.seq())?;
        }

    }
    Ok(())
}
