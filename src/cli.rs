use std::path::PathBuf;

use clap::{Parser, Subcommand};

#[derive(Parser)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,

    /// Input sequence file. Uncompressed or bgzipped.
    #[arg(short, long, global = true)]
    pub infile: Option<PathBuf>,

    /// Input bed file. Each region should map to a sequence from infile.
    #[arg(short = 'r', long, global = true)]
    pub inbedfile: Option<PathBuf>,

    /// Output sequence file.
    #[arg(short, long, global = true)]
    pub outfile: Option<PathBuf>,

    /// Output BED file with misassemblies.
    #[arg(short = 'b', long, global = true)]
    pub outbedfile: Option<PathBuf>,

    /// Seed to use for the random number generator.
    #[arg(short, long, global = true)]
    pub seed: Option<u64>,

    /// Randomize length.
    #[arg(long, action, default_value_t = false, global = true)]
    pub randomize_length: bool,
}

#[derive(Debug, PartialEq, Eq, Subcommand)]
pub enum Commands {
    /// Simulate a misjoin in a sequence.
    Misjoin {
        /// Number of misjoins to simulate.
        #[arg(short, long, default_value_t = 1)]
        number: usize,

        /// Max length of misjoin.
        #[arg(short, long, default_value_t = 5_000)]
        length: usize,
    },

    /// Simulate a falsely duplicated sequence.
    FalseDuplication {
        /// Number of false duplications to simulate.
        #[arg(short, long, default_value_t = 1)]
        number: usize,

        /// Max length of sequence to duplicate.
        #[arg(short, long, default_value_t = 5_000)]
        length: usize,

        /// Maximum number of duplications for any single segment.
        #[arg(short, long, default_value_t = 3)]
        max_duplications: usize,
    },

    /// Simulate a gap in a sequence.
    Gap {
        /// Number of gaps to simulate.
        #[arg(short, long, default_value_t = 1)]
        number: usize,

        /// Max length of gap simulate.
        #[arg(short, long, default_value_t = 5_000)]
        length: usize,
    },

    /// Simulate a break in a sequence.
    Break {
        /// Number of breaks to simulate.
        #[arg(short, long, default_value_t = 1)]
        number: usize,
    },
}
