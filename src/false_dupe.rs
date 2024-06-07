use eyre::ContextCompat;
use iset::IntervalSet;
use itertools::Itertools;
use noodles::{
    bed::{
        self,
        record::{Builder, OptionalFields},
    },
    core::Position,
};
use rand::{rngs::StdRng, seq::IteratorRandom, SeedableRng};

use crate::utils::generate_random_seq_ranges;

#[derive(Debug, PartialEq, Eq, Clone)]
pub struct DuplicateSequence {
    /// The duplicated sequence.
    pub seq: String,
    /// The duplicated segments.
    pub duplicated_seqs: Vec<Repeat>,
}

#[derive(Debug, PartialEq, Eq, Clone, Hash)]
pub struct Repeat {
    pub seq: String,
    pub start: usize,
    pub count: usize,
}

impl From<Repeat> for Builder<3> {
    fn from(rp: Repeat) -> Self {
        bed::Record::<3>::builder()
            .set_start_position(Position::new(rp.start.clamp(1, usize::MAX)).unwrap())
            .set_end_position(Position::new(rp.start + (rp.seq.len() * rp.count)).unwrap())
            .set_optional_fields(OptionalFields::from(vec![rp.count.to_string(), rp.seq]))
    }
}

pub fn generate_false_duplication(
    seq: &str,
    regions: &IntervalSet<Position>,
    length: usize,
    number: usize,
    max_duplications: usize,
    seed: Option<u64>,
) -> eyre::Result<DuplicateSequence> {
    let seq_segments = generate_random_seq_ranges(seq.len(), regions, length, number, seed)?
        .context("No sequence segments")?
        .collect_vec();
    let mut seq_iter = seq_segments.into_iter().peekable();
    let mut new_seq = String::new();
    let mut duplicated_seqs = vec![];

    // Add starting sequence before first position.
    if let Some((_, _, rrange)) = seq_iter.peek() {
        new_seq.push_str(&seq[..rrange.start]);
    };

    // TODO: Look into characteristics of false duplications. Probably not completely random.
    let mut rng = seed.map_or(StdRng::from_entropy(), StdRng::seed_from_u64);
    while let Some((_, _, rrange)) = seq_iter.next() {
        let num_dupes = (2..max_duplications.clamp(1, usize::MAX))
            .choose(&mut rng)
            .unwrap();
        let dup_seq = &seq[rrange.clone()];
        let repeat = Repeat {
            seq: dup_seq.to_string(),
            start: rrange.start,
            count: num_dupes,
        };

        for _ in 0..num_dupes {
            new_seq.push_str(dup_seq);
        }

        let remaining_seq = if let Some((_, _, next_rrange)) = seq_iter.peek() {
            &seq[rrange.end..next_rrange.start]
        } else {
            &seq[rrange.end..seq.len()]
        };
        new_seq.push_str(remaining_seq);
        duplicated_seqs.push(repeat);
    }

    Ok(DuplicateSequence {
        seq: new_seq,
        duplicated_seqs,
    })
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_generate_false_duplication() {
        let seq = "AAAGGCCCTTTTCCGGGGGAACTTCGGAC";
        let regions = IntervalSet::from_iter(std::iter::once(
            Position::new(1).unwrap()..Position::new(seq.len()).unwrap(),
        ));

        let new_seq = generate_false_duplication(seq, &regions, 10, 1, 3, Some(432)).unwrap();
        assert_eq!(
            new_seq,
            DuplicateSequence {
                seq: "AAAGGCCCTTTTCCGGGGGAACTTCGGATTCGGAC".to_string(),
                duplicated_seqs: [Repeat {
                    seq: "TTCGGA".to_string(),
                    start: 22,
                    count: 2
                }]
                .to_vec()
            }
        );
    }
}
