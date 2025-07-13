use rust_lapper::Interval;
use std::ops::Range;

use crate::sequence::{SequenceSegment, SequenceType};

pub fn create_breaks(
    regions: impl Iterator<Item = Interval<usize, Range<usize>>>,
) -> Vec<SequenceSegment> {
    regions
        .map(|itv| {
            SequenceSegment {
                typ: SequenceType::Break,
                itv: Interval {
                    start: itv.val.start,
                    // Needs to be null interval
                    stop: itv.val.end - 1,
                    val: None,
                },
            }
        })
        .collect()
}

#[cfg(test)]
mod test {
    use crate::sequence::SequenceType;
    use rust_lapper::Interval;

    use super::*;

    #[test]
    fn test_generate_breaks() {
        let seq = "AAAGGCCCGGCCCGGGGATTTTATTTTGGGCCGCCCAATTTAATTT";
        let regions = vec![
            Interval {
                start: 1,
                stop: seq.len(),
                val: 17..18,
            },
            Interval {
                start: 1,
                stop: seq.len(),
                val: 25..26,
            },
            Interval {
                start: 1,
                stop: seq.len(),
                val: 45..46,
            },
        ];

        let res = create_breaks(regions.into_iter());
        assert_eq!(
            res,
            [
                SequenceSegment {
                    typ: SequenceType::Break,
                    itv: Interval {
                        start: 17,
                        stop: 17,
                        val: None
                    }
                },
                SequenceSegment {
                    typ: SequenceType::Break,
                    itv: Interval {
                        start: 25,
                        stop: 25,
                        val: None
                    }
                },
                SequenceSegment {
                    typ: SequenceType::Break,
                    itv: Interval {
                        start: 45,
                        stop: 45,
                        val: None
                    }
                }
            ]
        )
    }
}
