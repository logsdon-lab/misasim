use eyre::ContextCompat;
use iset::IntervalSet;
use itertools::Itertools;
use noodles::{
    bed::{
        record::{Builder, OptionalFields},
        Record,
    },
    core::Position,
};

use crate::utils::generate_random_seq_ranges;

#[derive(Debug, PartialEq, Eq, Clone)]
pub struct RemovedSequence<'a> {
    pub start: usize,
    pub end: usize,
    pub seq: &'a str,
}

impl<'a> TryFrom<RemovedSequence<'a>> for Builder<3> {
    type Error = eyre::Error;

    fn try_from(rem_seq: RemovedSequence) -> Result<Self, eyre::Error> {
        Ok(Record::builder()
            .set_start_position(Position::new(rem_seq.start).context("Zero start position")?)
            .set_end_position(Position::new(rem_seq.end).context("Zero end position")?)
            .set_optional_fields(OptionalFields::from(vec![rem_seq.seq.to_owned()])))
    }
}

#[derive(Debug, PartialEq, Eq)]
pub struct DeletedSequence<'a> {
    pub seq: String,
    pub removed_seqs: Vec<RemovedSequence<'a>>,
}

pub fn generate_deletion<'a>(
    seq: &'a str,
    regions: &IntervalSet<Position>,
    length: usize,
    number_dels: usize,
    mask_del: bool,
    seed: Option<u64>,
    randomize_length: bool,
) -> eyre::Result<DeletedSequence<'a>> {
    let mut new_seq = String::with_capacity(seq.len());
    let mut removed_seqs: Vec<RemovedSequence> = Vec::with_capacity(number_dels);
    let seq_segments = generate_random_seq_ranges(
        seq.len(),
        regions,
        length,
        number_dels,
        seed,
        randomize_length,
    )?
    .context("No sequence segments")?
    .collect_vec();

    let mut seq_iter = seq_segments.into_iter().peekable();
    // Add starting sequence before first position.
    if let Some((_, _, del_range)) = seq_iter.peek() {
        new_seq.push_str(&seq[..del_range.start]);
    };

    while let Some((_, _, rrange)) = seq_iter.next() {
        let del_seq = &seq[rrange.clone()];
        if mask_del {
            new_seq.push_str(&"N".repeat(del_seq.len()));
        }

        removed_seqs.push(RemovedSequence {
            start: rrange.start,
            end: rrange.end,
            seq: del_seq,
        });

        let remaining_seq = if let Some((_, _, next_rrange)) = seq_iter.peek() {
            &seq[rrange.end..next_rrange.start]
        } else {
            &seq[rrange.end..seq.len()]
        };

        new_seq.push_str(remaining_seq);
    }

    Ok(DeletedSequence {
        seq: new_seq,
        removed_seqs,
    })
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_generate_misjoin() {
        let seq = "AAAGGCCCGGCCCGGGGATTTTATTTTGGGCCGCCCAATTTAATTT";
        let regions = IntervalSet::from_iter(std::iter::once(
            Position::new(1).unwrap()..Position::new(seq.len()).unwrap(),
        ));
        let new_seq = generate_deletion(seq, &regions, 10, 1, false, Some(42), true).unwrap();

        assert_eq!(
            DeletedSequence {
                seq: "AAAGGCCCGGCCCGGGGATTTTATGGGCCGCCCAATTTAATTT".to_string(),
                removed_seqs: [RemovedSequence {
                    start: 24,
                    end: 27,
                    seq: "TTT"
                }]
                .to_vec()
            },
            new_seq
        );
    }

    #[test]
    fn test_generate_misjoin_multiple() {
        let seq = "AAAGGCCCGGCCCGGGGATTTTATTTTGGGCCGCCCAATTTAATTT";
        let regions = IntervalSet::from_iter(std::iter::once(
            Position::new(1).unwrap()..Position::new(seq.len()).unwrap(),
        ));
        let new_seq = generate_deletion(seq, &regions, 10, 3, false, Some(42), true).unwrap();

        assert_eq!(
            DeletedSequence {
                seq: "AAAGGCCCGGCCCGGGGGGCCGCCCAATTTAATT".to_string(),
                removed_seqs: [
                    RemovedSequence {
                        start: 16,
                        end: 24,
                        seq: "GATTTTAT"
                    },
                    RemovedSequence {
                        start: 24,
                        end: 27,
                        seq: "TTT"
                    },
                    RemovedSequence {
                        start: 44,
                        end: 45,
                        seq: "T"
                    }
                ]
                .to_vec()
            },
            new_seq
        );
    }

    #[test]
    fn test_generate_gap_multiple() {
        let seq = "AAAGGCCCGGCCCGGGGATTTTATTTTGGGCCGCCCAATTTAATTT";
        let regions = IntervalSet::from_iter(std::iter::once(
            Position::new(1).unwrap()..Position::new(seq.len()).unwrap(),
        ));
        let new_seq = generate_deletion(seq, &regions, 10, 3, true, Some(42), true).unwrap();

        assert_eq!(
            DeletedSequence {
                seq: "AAAGGCCCGGCCCGGGNNNNNNNNNNNGGGCCGCCCAATTTAATNT".to_string(),
                removed_seqs: [
                    RemovedSequence {
                        start: 16,
                        end: 24,
                        seq: "GATTTTAT"
                    },
                    RemovedSequence {
                        start: 24,
                        end: 27,
                        seq: "TTT"
                    },
                    RemovedSequence {
                        start: 44,
                        end: 45,
                        seq: "T"
                    }
                ]
                .to_vec()
            },
            new_seq
        )
    }
}
