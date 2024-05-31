use std::path::PathBuf;

use clap::{Parser, Subcommand};

#[derive(Parser)]
pub struct Cli {
    /// Input sequence file.
    pub infile: PathBuf,

    /// Output sequence file.
    pub outfile: Option<PathBuf>,

    /// Output BED file with misassemblies.
    pub outbedfile: Option<PathBuf>,

    /// Seed to use for the random number generator.
    #[arg(short, long)]
    pub seed: Option<u64>,

    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Simulate a misjoin in a sequence.
    Misjoin {
        /// Length in bp to remove.
        #[arg(short, long)]
        length: usize,

        /// Number of misjoins to simulate.
        #[arg(short, long)]
        number: usize,
    },

    /// Simulate a collapse in a sequence.
    Collapse {
        /// Length in bp to collapse.
        /// Requires a repetitive sequence of this length or more in the input.
        #[arg(short, long)]
        length: usize,

        /// Number of repeats to collapse.
        #[arg(short, long)]
        num_repeats: usize,
    },

    /// Simulate a falsely duplicated sequence.
    FalseDuplication {
        /// Length in bp to duplicate.
        #[arg(short, long)]
        length: usize,

        /// Number of false duplications to simulate.
        #[arg(short, long)]
        number: usize,
    },

    /// Simulate a gap in a sequence.
    Gap {
        /// Length in bp to insert.
        #[arg(short, long)]
        length: usize,

        /// Number of gaps to simulate.
        #[arg(short, long)]
        number: usize,
    }
}