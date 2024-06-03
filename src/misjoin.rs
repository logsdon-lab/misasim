use eyre::ContextCompat;
use noodles::{
    bed::{
        record::{Builder, OptionalFields},
        Record,
    },
    core::Position,
};
use rand::prelude::*;

#[derive(Debug, PartialEq, Eq, Clone)]
pub struct DeletedSequence<'a> {
    pub start: usize,
    pub end: usize,
    pub seq: &'a str,
}

impl<'a> TryFrom<DeletedSequence<'a>> for Builder<3> {
    type Error = eyre::Error;

    fn try_from(del_seq: DeletedSequence) -> Result<Self, eyre::Error> {
        Ok(Record::builder()
            .set_start_position(Position::new(del_seq.start).context("Zero start position")?)
            .set_end_position(Position::new(del_seq.end).context("Zero end position")?)
            .set_optional_fields(OptionalFields::from(vec![del_seq.seq.to_owned()])))
    }
}

#[derive(Debug, PartialEq, Eq)]
pub struct MisjoinedSequence<'a> {
    pub seq: String,
    pub del_seqs: Vec<DeletedSequence<'a>>,
}

// TODO: Add length. Difficult as may overlap with another. Ignoring for now.
pub fn generate_misjoin(
    seq: &str,
    num_misjoins: usize,
    seed: Option<u64>,
) -> eyre::Result<MisjoinedSequence> {
    let mut new_seq = String::with_capacity(seq.len());
    let mut del_seqs: Vec<DeletedSequence> = Vec::with_capacity(num_misjoins);
    let mut rng = seed.map_or(StdRng::from_entropy(), StdRng::seed_from_u64);
    let mut misjoin_pos = (0..seq.len()).choose_multiple(&mut rng, num_misjoins);
    misjoin_pos.sort();

    let bp_dels = misjoin_pos.iter().enumerate().map(|(i, p)| {
        let (stop_pos, dst) = if let Some(pos) = misjoin_pos.get(i + 1) {
            (*pos, pos - p)
        } else {
            (seq.len(), seq.len() - p)
        };
        (stop_pos, rng.gen_range(0..dst))
    });

    // Add starting sequence before first position.
    if let Some(pos) = misjoin_pos.first() {
        new_seq.push_str(&seq[..*pos]);
    };

    for (pos, (stop, bp_del)) in misjoin_pos.iter().zip(bp_dels) {
        let remaining_seq = &seq[pos + bp_del..stop];

        let del_start = *pos;
        let del_end = *pos + bp_del;
        del_seqs.push(DeletedSequence {
            start: del_start,
            end: del_end,
            seq: &seq[del_start..del_end],
        });
        new_seq.push_str(remaining_seq);
    }

    Ok(MisjoinedSequence {
        seq: new_seq,
        del_seqs,
    })
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_generate_misjoin() {
        let seq = "AAAGGCCCGGCCCGGGGATTTTATTTTGGGCCGCCCAATTTAATTT";
        let new_seq = generate_misjoin(seq, 1, Some(42)).unwrap();

        assert_eq!(
            MisjoinedSequence {
                seq: "AAAGGCCCGGCCCGGGGATTTTAATTT".to_string(),
                del_seqs: [DeletedSequence {
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
        let new_seq = generate_misjoin(seq, 3, Some(42)).unwrap();

        assert_eq!(
            MisjoinedSequence {
                seq: "AAAGGGGATTTTATTTTGGCT".to_string(),
                del_seqs: [
                    DeletedSequence {
                        start: 4,
                        end: 14,
                        seq: "GCCCGGCCCG"
                    },
                    DeletedSequence {
                        start: 29,
                        end: 34,
                        seq: "GCCGC"
                    },
                    DeletedSequence {
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
}
