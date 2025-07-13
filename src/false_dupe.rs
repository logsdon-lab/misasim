use std::ops::Range;

use rand::{rngs::StdRng, seq::IteratorRandom, SeedableRng};
use rust_lapper::Interval;

use crate::sequence::{SequenceSegment, SequenceType};

pub fn create_false_dupe(
    seq: &str,
    regions: impl Iterator<Item = Interval<usize, Range<usize>>>,
    seed: Option<u64>,
    max_duplications: usize,
) -> Vec<SequenceSegment> {
    let mut rng = seed.map_or(StdRng::from_entropy(), StdRng::seed_from_u64);
    regions
        .into_iter()
        .map(
            move |Interval {
                      start: _,
                      stop: _,
                      val: range,
                  }| {
                let num_dupes = (2..(max_duplications + 1).clamp(3, usize::MAX))
                    .choose(&mut rng)
                    .expect("Require at minimum 2 for max duplications.");
                let dup_seq = &seq[range.clone()];
                SequenceSegment {
                    typ: SequenceType::FalseDuplication,
                    itv: Interval {
                        start: range.start,
                        stop: range.end,
                        val: Some(dup_seq.repeat(num_dupes)),
                    },
                }
            },
        )
        .collect()
}

#[cfg(test)]
mod test {
    use rust_lapper::Interval;

    use super::*;

    #[test]
    fn test_generate_false_duplication() {
        let seq = "AAAGGCCCTTTTCCGGGGGAACTTCGGAC";
        let regions = vec![Interval {
            start: 1,
            stop: seq.len(),
            val: 22..28,
        }];

        let res = create_false_dupe(seq, regions.into_iter(), Some(10), 2);
        assert_eq!(
            res,
            vec![SequenceSegment {
                typ: SequenceType::FalseDuplication,
                itv: Interval {
                    start: 22,
                    stop: 28,
                    val: Some("TTCGGATTCGGA".to_owned(),)
                }
            }]
        );
    }
}
