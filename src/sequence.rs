use rust_lapper::Interval;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SequenceType {
    Good,
    Misjoin,
    FalseDuplication,
    Gap,
    Break,
    Inversion,
}

impl SequenceType {
    pub fn as_color(&self) -> &'static str {
        match self {
            SequenceType::Good => "0,0,0",
            SequenceType::Misjoin
            | SequenceType::FalseDuplication
            | SequenceType::Gap
            | SequenceType::Break
            | SequenceType::Inversion => "255,0,0",
        }
    }
}

impl std::fmt::Display for SequenceType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let seq_type = match self {
            SequenceType::Good => "good",
            SequenceType::Misjoin => "misjoin",
            SequenceType::FalseDuplication => "false_dupe",
            SequenceType::Gap => "gap",
            SequenceType::Break => "break",
            SequenceType::Inversion => "inversion",
        };
        write!(f, "{seq_type}")
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SequenceSegment {
    pub typ: SequenceType,
    pub itv: Interval<usize, Option<String>>,
}
