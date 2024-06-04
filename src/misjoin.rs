use eyre::ContextCompat;
use itertools::Itertools;
use noodles::{
    bed::{
        record::{Builder, OptionalFields},
        Record,
    },
    core::Position,
};

use crate::utils::get_sequence_segments;

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

// TODO: Add length. Difficult as may overlap with another. Ignoring for now.
pub fn generate_deletion(
    seq: &str,
    number_dels: usize,
    mask_del: bool,
    seed: Option<u64>,
) -> eyre::Result<DeletedSequence> {
    let mut new_seq = String::with_capacity(seq.len());
    let mut removed_seqs: Vec<RemovedSequence> = Vec::with_capacity(number_dels);
    let seq_segments = get_sequence_segments(seq, number_dels, seed)
        .context("No sequence segments")?
        .collect_vec();

    // Add starting sequence before first position.
    if let Some((start, _, _)) = seq_segments.first() {
        new_seq.push_str(&seq[..*start]);
    };

    for (pos, _, rrange) in seq_segments {
        if mask_del {
            new_seq.push_str(&"N".repeat(rrange.len()));
        }
        let del_start = pos;
        let del_end = pos + rrange.len();
        let del_seq = &seq[del_start..del_end];
        let remaining_seq = &seq[rrange];
        new_seq.push_str(remaining_seq);

        removed_seqs.push(RemovedSequence {
            start: del_start,
            end: del_end,
            seq: del_seq,
        });
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
        let new_seq = generate_deletion(seq, 1, false, Some(42)).unwrap();

        assert_eq!(
            DeletedSequence {
                seq: "AAAGGCCCGGCCCGGGGATTTTAATTT".to_string(),
                removed_seqs: [RemovedSequence {
                    start: 22,
                    end: 41,
                    seq: "ATTTTGGGCCGCCCAATTT"
                }]
                .to_vec()
            },
            new_seq
        );
    }

    #[test]
    fn test_generate_misjoin_multiple() {
        let seq = "AAAGGCCCGGCCCGGGGATTTTATTTTGGGCCGCCCAATTTAATTT";
        let new_seq = generate_deletion(seq, 3, false, Some(42)).unwrap();

        assert_eq!(
            DeletedSequence {
                seq: "AAAGGGGATTTTATTTTGGCT".to_string(),
                removed_seqs: [
                    RemovedSequence {
                        start: 4,
                        end: 14,
                        seq: "GCCCGGCCCG"
                    },
                    RemovedSequence {
                        start: 29,
                        end: 34,
                        seq: "GCCGC"
                    },
                    RemovedSequence {
                        start: 35,
                        end: 45,
                        seq: "CAATTTAATT"
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
        let new_seq = generate_deletion(seq, 3, true, Some(42)).unwrap();

        assert_eq!(
            DeletedSequence {
                seq: "AAAGNNNNNNNNNNGGGATTTTATTTTGGNNNNNCNNNNNNNNNNT".to_string(),
                removed_seqs: [
                    RemovedSequence {
                        start: 4,
                        end: 14,
                        seq: "GCCCGGCCCG"
                    },
                    RemovedSequence {
                        start: 29,
                        end: 34,
                        seq: "GCCGC"
                    },
                    RemovedSequence {
                        start: 35,
                        end: 45,
                        seq: "CAATTTAATT"
                    }
                ]
                .to_vec()
            },
            new_seq
        )
    }
}
