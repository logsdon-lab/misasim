use std::ops::Range;

use rust_lapper::Interval;

use crate::sequence::{SequenceSegment, SequenceType};

pub fn create_deletion(
    seq: &str,
    regions: impl Iterator<Item = Interval<usize, Range<usize>>>,
    mask: bool,
) -> Vec<SequenceSegment> {
    regions
        .into_iter()
        .map(
            move |Interval {
                      start: _,
                      stop: _,
                      val: range,
                  }| {
                let del_seq = &seq[range.clone()];
                let (seq_typ, del_seq) = if mask {
                    (SequenceType::Gap, Some("N".repeat(del_seq.len())))
                } else {
                    (SequenceType::Misjoin, None)
                };
                SequenceSegment {
                    typ: seq_typ,
                    itv: Interval {
                        start: range.start,
                        stop: range.end,
                        val: del_seq,
                    },
                }
            },
        )
        .collect()
}

#[cfg(test)]
mod test {
    use crate::sequence::{SequenceSegment, SequenceType};
    use rust_lapper::Interval;

    use super::*;

    #[test]
    fn test_generate_misjoin() {
        let seq = "AAAGGCCCGGCCCGGGGATTTTATTTTGGGCCGCCCAATTTAATTT";
        let regions = vec![Interval {
            start: 1,
            stop: seq.len(),
            val: 24..27,
        }];

        let res = create_deletion(seq, regions.into_iter(), false);

        assert_eq!(
            res,
            vec![SequenceSegment {
                typ: SequenceType::Misjoin,
                itv: Interval {
                    start: 24,
                    stop: 27,
                    val: None
                }
            }]
        );
    }

    #[test]
    fn test_generate_gap_multiple() {
        let seq = "AAAGGCCCGGCCCGGGGATTTTATTTTGGGCCGCCCAATTTAATTT";
        let regions = vec![
            Interval {
                start: 1,
                stop: seq.len(),
                val: 16..24,
            },
            Interval {
                start: 1,
                stop: seq.len(),
                val: 24..27,
            },
            Interval {
                start: 1,
                stop: seq.len(),
                val: 44..45,
            },
        ];
        let res = create_deletion(seq, regions.into_iter(), true);

        // AAAGGCCCGGCCCGGGNNNNNNNNNNNGGGCCGCCCAATTTAATNT
        assert_eq!(
            res,
            vec![
                SequenceSegment {
                    typ: SequenceType::Gap,
                    itv: Interval {
                        start: 16,
                        stop: 24,
                        val: Some("NNNNNNNN".to_owned())
                    }
                },
                SequenceSegment {
                    typ: SequenceType::Gap,
                    itv: Interval {
                        start: 24,
                        stop: 27,
                        val: Some("NNN".to_owned())
                    }
                },
                SequenceSegment {
                    typ: SequenceType::Gap,
                    itv: Interval {
                        start: 44,
                        stop: 45,
                        val: Some("N".to_owned())
                    }
                },
            ]
        )
    }
}
