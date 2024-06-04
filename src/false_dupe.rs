use eyre::ContextCompat;
use itertools::Itertools;
use rand::{rngs::StdRng, seq::IteratorRandom, SeedableRng};

use crate::{collapse::Repeat, utils::get_sequence_segments};

#[derive(Debug, PartialEq, Eq, Clone)]
pub struct DuplicateSequence {
    /// The duplicated sequence.
    pub seq: String,
    /// The duplicated segments.
    pub duplicated_seqs: Vec<Repeat>,
}

pub fn generate_false_duplication(
    seq: &str,
    number: usize,
    max_duplications: usize,
    seed: Option<u64>,
) -> eyre::Result<DuplicateSequence> {
    let seq_segments = get_sequence_segments(seq, number, seed)
        .context("No sequence segments")?
        .collect_vec();
    let mut new_seq = String::new();
    let mut duplicated_seqs = vec![];

    // Add starting sequence before first position.
    if let Some((start, _, _)) = seq_segments.first() {
        new_seq.push_str(&seq[..*start]);
    };

    // TODO: Look into characteristics of false duplications. Probably not completely random.
    let mut rng = seed.map_or(StdRng::from_entropy(), StdRng::seed_from_u64);
    for (start, _, rrange) in seq_segments.into_iter() {
        let num_dupes = (1..max_duplications.clamp(1, usize::MAX))
            .choose(&mut rng)
            .unwrap();
        let dup_seq = &seq[start..start + rrange.len()];
        let repeat = Repeat {
            seq: dup_seq.to_string(),
            start,
            count: num_dupes,
        };

        for _ in 0..num_dupes {
            new_seq.push_str(dup_seq);
        }
        new_seq.push_str(&seq[rrange]);
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
        let new_seq = generate_false_duplication(seq, 1, 3, Some(432)).unwrap();
        assert_eq!(
            new_seq,
            DuplicateSequence {
                seq: "AAAGGCCCTTTTCCGGGGGAACTTCGGACGGGAACTTCGGACGGGAACTTCGGAC".to_string(),
                duplicated_seqs: [Repeat {
                    seq: "GGGAACTTCGGAC".to_string(),
                    start: 16,
                    count: 2
                }]
                .to_vec()
            }
        );
    }
}
