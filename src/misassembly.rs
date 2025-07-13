use std::{fmt::Debug, fs::read_to_string};

use eyre::bail;
use json::JsonValue;
use rust_lapper::{Interval, Lapper};

use crate::{
    breaks::create_breaks,
    cli::Misassembly,
    false_dupe::create_false_dupe,
    inversion::create_inversion,
    misjoin::create_deletion,
    sequence::SequenceSegment,
    utils::{generate_random_seq_ranges, subtract_intervals, subtract_misassembled_sequences},
};

impl Misassembly {
    /// Generate split sequences from `seq` within given regions.
    fn generate_split_sequences<T>(
        &self,
        seq: &str,
        regions: &Lapper<usize, T>,
        seed: Option<u64>,
        randomize_length: bool,
    ) -> eyre::Result<Vec<SequenceSegment>>
    where
        T: Eq + Clone + Send + Sync + Debug,
    {
        Ok(match self {
            Misassembly::Misjoin { number, length } => create_deletion(
                seq,
                generate_random_seq_ranges(
                    seq.len(),
                    regions,
                    *length,
                    *number,
                    seed,
                    randomize_length,
                )?,
                false,
            ),
            Misassembly::Gap { number, length } => create_deletion(
                seq,
                generate_random_seq_ranges(
                    seq.len(),
                    regions,
                    *length,
                    *number,
                    seed,
                    randomize_length,
                )?,
                true,
            ),
            Misassembly::Inversion { number, length } => create_inversion(
                seq,
                generate_random_seq_ranges(
                    seq.len(),
                    regions,
                    *length,
                    *number,
                    seed,
                    randomize_length,
                )?,
            ),
            Misassembly::FalseDuplication {
                number,
                length,
                max_duplications,
            } => create_false_dupe(
                seq,
                generate_random_seq_ranges(
                    seq.len(),
                    regions,
                    *length,
                    *number,
                    seed,
                    randomize_length,
                )?,
                seed,
                *max_duplications,
            ),
            Misassembly::Break { number } => create_breaks(generate_random_seq_ranges(
                seq.len(),
                regions,
                1,
                *number,
                seed,
                false,
            )?),
            // TODO: Haplotype switch.
            _ => bail!("Invalid option. {self:?}"),
        })
    }
}

pub fn create_all_sequences(
    command: &Misassembly,
    seq: &str,
    record_regions: &Lapper<usize, ()>,
    seed: Option<u64>,
    randomize_length: bool,
) -> eyre::Result<Vec<SequenceSegment>> {
    let mut all_sequences = vec![];

    if let Misassembly::Multiple { path } = command {
        let mut original_record_regions = record_regions.clone();

        let json_str = read_to_string(path)?;
        let JsonValue::Array(json_misassemblies) = json::parse(&json_str)? else {
            bail!("First JSON value must be an array. ({json_str})")
        };

        for (i, json_misassembly) in json_misassemblies.into_iter().enumerate() {
            // Offset seed by 1 to ensure locations are different.
            let seed = seed.map(|seed| seed + i as u64);

            let misassembly: Misassembly = json_misassembly.try_into()?;
            let misassembly_seqs = misassembly.generate_split_sequences(
                seq,
                &original_record_regions,
                seed,
                randomize_length,
            )?;
            all_sequences.extend(misassembly_seqs.iter().cloned());
            let itree_seqs = Lapper::new(misassembly_seqs.into_iter().map(|s| s.itv).collect());

            // Then remove these regions from lapper.
            let mut new_itvs = vec![];
            for itv in original_record_regions.iter() {
                let ovl = itree_seqs.find(itv.start, itv.stop);
                new_itvs.extend(subtract_intervals(
                    itv.clone(),
                    ovl.into_iter().map(|itv| Interval {
                        start: itv.start,
                        stop: itv.stop,
                        val: (),
                    }),
                ));
            }
            original_record_regions = Lapper::new(new_itvs);
        }
    } else {
        let mut seqs =
            command.generate_split_sequences(seq, record_regions, seed, randomize_length)?;
        all_sequences.append(&mut seqs);
    }

    let good_sequence_itvs = subtract_misassembled_sequences(seq, all_sequences.iter());
    all_sequences.extend(good_sequence_itvs);

    all_sequences.sort_by(|a, b| a.itv.start.cmp(&b.itv.start));

    Ok(all_sequences)
}

#[cfg(test)]
mod test {
    use crate::{
        cli::Misassembly,
        misassembly::create_all_sequences,
        sequence::{SequenceSegment, SequenceType},
    };
    use rust_lapper::{Interval, Lapper};
    use std::{path::PathBuf, str::FromStr};

    #[test]
    fn test_all_false_duplication() {
        let seq = "ATTATTATTGCA";
        let res = create_all_sequences(
            &Misassembly::FalseDuplication {
                number: 1,
                length: 2,
                max_duplications: 2,
            },
            seq,
            &Lapper::new(vec![Interval {
                start: 1,
                stop: 5,
                val: (),
            }]),
            Some(12),
            true,
        )
        .unwrap();

        assert_eq!(
            res,
            vec![
                SequenceSegment {
                    typ: SequenceType::Good,
                    itv: Interval {
                        start: 0,
                        stop: 3,
                        val: None
                    }
                },
                SequenceSegment {
                    typ: SequenceType::FalseDuplication,
                    itv: Interval {
                        start: 3,
                        stop: 5,
                        val: Some("ATAT".to_owned())
                    }
                },
                SequenceSegment {
                    typ: SequenceType::Good,
                    itv: Interval {
                        start: 5,
                        stop: 12,
                        val: None
                    }
                },
            ]
        )
    }

    #[test]
    fn test_all_breaks() {
        let seq = "ATTATTATTGCA";
        let res = create_all_sequences(
            &Misassembly::Break { number: 2 },
            seq,
            &Lapper::new(vec![Interval {
                start: 1,
                stop: seq.len(),
                val: (),
            }]),
            Some(12),
            true,
        )
        .unwrap();

        assert_eq!(
            res,
            [
                SequenceSegment {
                    typ: SequenceType::Good,
                    itv: Interval {
                        start: 0,
                        stop: 2,
                        val: None
                    }
                },
                SequenceSegment {
                    typ: SequenceType::Break,
                    itv: Interval {
                        start: 2,
                        stop: 2,
                        val: None
                    }
                },
                SequenceSegment {
                    typ: SequenceType::Good,
                    itv: Interval {
                        start: 2,
                        stop: 10,
                        val: None
                    }
                },
                SequenceSegment {
                    typ: SequenceType::Break,
                    itv: Interval {
                        start: 10,
                        stop: 10,
                        val: None
                    }
                },
                SequenceSegment {
                    typ: SequenceType::Good,
                    itv: Interval {
                        start: 10,
                        stop: 12,
                        val: None
                    }
                }
            ]
        )
    }

    #[test]
    fn test_all_misjoin() {
        let seq = "ATTATTATTGCA";
        let res = create_all_sequences(
            &Misassembly::Misjoin {
                number: 2,
                length: 4,
            },
            seq,
            &Lapper::new(vec![Interval {
                start: 1,
                stop: 5,
                val: (),
            }]),
            Some(12),
            true,
        )
        .unwrap();

        assert_eq!(
            res,
            vec![
                SequenceSegment {
                    typ: SequenceType::Good,
                    itv: Interval {
                        start: 0,
                        stop: 2,
                        val: None
                    }
                },
                SequenceSegment {
                    typ: SequenceType::Misjoin,
                    itv: Interval {
                        start: 2,
                        stop: 3,
                        val: None
                    }
                },
                SequenceSegment {
                    typ: SequenceType::Misjoin,
                    itv: Interval {
                        start: 3,
                        stop: 6,
                        val: None
                    }
                },
                SequenceSegment {
                    typ: SequenceType::Good,
                    itv: Interval {
                        start: 6,
                        stop: 12,
                        val: None
                    }
                }
            ]
        )
    }

    #[test]
    fn test_all_multiple_types() {
        let seq = "ATTATTATTGCA";
        let file = "test/data/multiple.json";
        let res = create_all_sequences(
            &Misassembly::Multiple {
                path: PathBuf::from_str(file).unwrap(),
            },
            seq,
            &Lapper::new(vec![Interval {
                start: 1,
                stop: seq.len(),
                val: (),
            }]),
            Some(12),
            true,
        )
        .unwrap();

        assert_eq!(
            res,
            vec![
                SequenceSegment {
                    typ: SequenceType::Good,
                    itv: Interval {
                        start: 0,
                        stop: 2,
                        val: None
                    }
                },
                SequenceSegment {
                    typ: SequenceType::Inversion,
                    itv: Interval {
                        start: 2,
                        stop: 4,
                        val: Some("TA".to_owned())
                    }
                },
                SequenceSegment {
                    typ: SequenceType::Good,
                    itv: Interval {
                        start: 4,
                        stop: 6,
                        val: None
                    }
                },
                SequenceSegment {
                    typ: SequenceType::Misjoin,
                    itv: Interval {
                        start: 6,
                        stop: 9,
                        val: None
                    }
                },
                SequenceSegment {
                    typ: SequenceType::Good,
                    itv: Interval {
                        start: 9,
                        stop: 12,
                        val: None
                    }
                },
            ]
        )
    }
}
