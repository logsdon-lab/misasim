use rust_lapper::Interval;
use std::ops::Range;

use crate::sequence::{SequenceSegment, SequenceType};

pub fn create_inversion(
    seq: &str,
    regions: impl Iterator<Item = Interval<usize, Range<usize>>>,
) -> Vec<SequenceSegment> {
    regions
        .into_iter()
        .map(
            |Interval {
                 start: _,
                 stop: _,
                 val: range,
             }| {
                let inv_seq: String = seq[range.clone()]
                    .chars()
                    .map(|nt| match nt.to_ascii_uppercase() {
                        'A' => 'T',
                        'T' => 'A',
                        'G' => 'C',
                        'C' => 'G',
                        _ => nt,
                    })
                    .rev()
                    .collect();

                SequenceSegment {
                    typ: SequenceType::Inversion,
                    itv: Interval {
                        start: range.start,
                        stop: range.end,
                        val: Some(inv_seq),
                    },
                }
            },
        )
        .collect()
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_generate_inversion() {
        let seq = "AAAGGCCCGGCCCGGGGATTTTATTTTGGGCCGCCCAATTTAATTT";
        //                      "--------------------------------------TTTAAT--";
        let regions = vec![Interval {
            start: 1,
            stop: seq.len(),
            val: 38..44,
        }];

        let res = create_inversion(seq, regions.into_iter());
        assert_eq!(
            res,
            vec![SequenceSegment {
                typ: SequenceType::Inversion,
                itv: Interval {
                    start: 38,
                    stop: 44,
                    val: Some("ATTAAA".to_owned())
                }
            }]
        )
    }
}
