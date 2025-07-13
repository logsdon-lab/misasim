use std::{collections::HashMap, path::PathBuf};

use clap::{Parser, Subcommand};
use eyre::bail;
use json::JsonValue;

use crate::sequence::SequenceType;

const DEFAULT_NUMBER: usize = 1;
const DEFAULT_LENGTH: usize = 5000;
const DEFAULT_FALSE_DUPE_MAX: usize = 2;

#[derive(Parser)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Misassembly,

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

    /// Group by regex pattern.
    /// ex. "^.*?_(?<hap>.*?)$" with group by haplotype.
    #[arg(short, long, global = true)]
    pub group_by: Option<String>,
}

#[derive(Debug, PartialEq, Eq, Subcommand)]
pub enum Misassembly {
    /// Simulate a misjoin in a sequence.
    Misjoin {
        /// Number of misjoins to simulate.
        #[arg(short, long, default_value_t = DEFAULT_NUMBER)]
        number: usize,

        /// Max length of misjoin.
        #[arg(short, long, default_value_t = DEFAULT_LENGTH)]
        length: usize,
    },

    /// Simulate a falsely duplicated sequence.
    FalseDuplication {
        /// Number of false duplications to simulate.
        #[arg(short, long, default_value_t = DEFAULT_NUMBER)]
        number: usize,

        /// Max length of sequence to duplicate.
        #[arg(short, long, default_value_t = DEFAULT_LENGTH)]
        length: usize,

        /// Maximum number of duplications for any single segment.
        #[arg(short, long, default_value_t = DEFAULT_FALSE_DUPE_MAX)]
        max_duplications: usize,
    },

    /// Simulate a gap in a sequence.
    Gap {
        /// Number of gaps to simulate.
        #[arg(short, long, default_value_t = DEFAULT_NUMBER)]
        number: usize,

        /// Max length of gap simulate.
        #[arg(short, long, default_value_t = DEFAULT_LENGTH)]
        length: usize,
    },

    /// Simulate a break in a sequence.
    Break {
        /// Number of breaks to simulate.
        #[arg(short, long, default_value_t = DEFAULT_NUMBER)]
        number: usize,
    },

    /// Simulate an inversion in a sequence.
    Inversion {
        /// Number of inversions to simulate.
        #[arg(short, long, default_value_t = DEFAULT_NUMBER)]
        number: usize,

        /// Max length of inversion simulate.
        #[arg(short, long, default_value_t = DEFAULT_LENGTH)]
        length: usize,
    },

    /// Simulate multiple misassembly types from an input JSON file.
    ///
    /// ex. A JSON file with a break and inversion.
    /// ```json
    /// [
    ///     {
    ///         "type": "break",
    ///         "number": 2
    ///     },
    ///     {
    ///         "type": "inversion",
    ///         "length": 5000
    ///     }
    /// ]
    /// ```
    Multiple {
        /// Path to JSON file.
        #[arg(short, long, required = true)]
        path: PathBuf,
    },
}

impl From<&Misassembly> for SequenceType {
    fn from(value: &Misassembly) -> Self {
        match value {
            Misassembly::Misjoin { .. } => Self::Misjoin,
            Misassembly::FalseDuplication { .. } => Self::FalseDuplication,
            Misassembly::Gap { .. } => Self::Gap,
            Misassembly::Break { .. } => Self::Break,
            Misassembly::Inversion { .. } => Self::Inversion,
            _ => Self::Good,
        }
    }
}

impl TryFrom<JsonValue> for Misassembly {
    type Error = eyre::Report;

    fn try_from(value: JsonValue) -> Result<Self, Self::Error> {
        let values: HashMap<&str, &JsonValue> = value.entries().collect();

        let Some(mtype) = values.get("mtype").and_then(|v| v.as_str()) else {
            bail!("Misassembly type not provided for {values:?}")
        };

        let number = values
            .get("number")
            .and_then(|v| v.as_usize())
            .unwrap_or(DEFAULT_NUMBER);
        let length = values
            .get("length")
            .and_then(|v| v.as_usize())
            .unwrap_or(DEFAULT_LENGTH);

        Ok(match mtype.to_lowercase().as_str() {
            "misjoin" => Misassembly::Misjoin { number, length },
            "false_duplication" => Misassembly::FalseDuplication {
                number,
                length,
                max_duplications: values
                    .get("max_duplications")
                    .and_then(|v| v.as_usize())
                    .unwrap_or(DEFAULT_FALSE_DUPE_MAX),
            },
            "gap" => Misassembly::Gap { number, length },
            "break" => Misassembly::Break { number },
            "inversion" => Misassembly::Inversion { number, length },
            _ => {
                bail!("Invalid mtype {mtype}")
            }
        })
    }
}

#[cfg(test)]
mod test {
    use json::object;

    use crate::cli::Misassembly;

    #[test]
    fn test_json_to_misassembly() {
        let cfg = object! {
            mtype: "misjoin",
            number: 1,
            length: 10
        };
        let res: Misassembly = cfg.try_into().unwrap();
        assert_eq!(
            res,
            Misassembly::Misjoin {
                number: 1,
                length: 10
            }
        )
    }

    #[test]
    fn test_json_to_misassembly_defaults() {
        let cfg = object! {
            mtype: "misjoin",
        };
        let res: Misassembly = cfg.try_into().unwrap();
        assert_eq!(
            res,
            Misassembly::Misjoin {
                number: 1,
                length: 5000
            }
        )
    }

    #[test]
    #[should_panic]
    fn test_json_to_misassembly_missing_mtype() {
        let cfg = object! {
            number: 1,
        };
        let _res: Misassembly = cfg.try_into().unwrap();
    }
}
