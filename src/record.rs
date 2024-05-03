use crate::{Error, ErrorKind, Result};
use std::{
    fmt::{Display, Formatter},
    str::FromStr,
};

/// A record in a HMMER tblout file
pub struct Record {
    target_name: String,
    query_name: String,
    hmm_from: i32,
    hmm_to: i32,
    ali_from: i32,
    ali_to: i32,
    env_from: i32,
    env_to: i32,
    sq_len: i32,
    strand: Strand,
    e_value: f32,
    score: f32,
    bias: f32,
}

impl Record {
    pub fn new(
        target_name: String,
        query_name: String,
        hmm_from: i32,
        hmm_to: i32,
        ali_from: i32,
        ali_to: i32,
        env_from: i32,
        env_to: i32,
        sq_len: i32,
        strand: Strand,
        e_value: f32,
        score: f32,
        bias: f32,
    ) -> Self {
        Record {
            target_name,
            query_name,
            hmm_from,
            hmm_to,
            ali_from,
            ali_to,
            env_from,
            env_to,
            sq_len,
            strand,
            e_value,
            score,
            bias,
        }
    }

    pub fn target_name(&self) -> String {
        self.target_name.clone()
    }

    pub fn query_name(&self) -> String {
        self.query_name.clone()
    }

    pub fn hmm_from(&self) -> i32 {
        self.hmm_from
    }
    pub fn hmm_to(&self) -> i32 {
        self.hmm_to
    }

    pub fn ali_from(&self) -> i32 {
        self.ali_from
    }

    pub fn ali_to(&self) -> i32 {
        self.ali_to
    }

    pub fn env_from(&self) -> i32 {
        self.env_from
    }

    pub fn env_to(&self) -> i32 {
        self.env_to
    }

    pub fn sq_len(&self) -> i32 {
        self.sq_len
    }

    pub fn strand(&self) -> Strand {
        self.strand
    }

    pub fn e_value(&self) -> f32 {
        self.e_value
    }

    pub fn score(&self) -> f32 {
        self.score
    }

    pub fn bias(&self) -> f32 {
        self.bias
    }
}

/// The strandedness of the HMM hit in the genome.
#[derive(Debug, Hash, PartialEq, Eq, Clone, Copy)]
pub enum Strand {
    Positive,
    Negative,
}

/// Implementation of `FromStr` for `Strand`.
impl FromStr for Strand {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        match s {
            "+" => Ok(Strand::Positive),
            "-" => Ok(Strand::Negative),
            _ => Err(Error::new(ErrorKind::ReadRecord(
                "The input was neither `-` nor `+`.".into(),
            ))),
        }
    }
}

/// An implementation of `Display` for `Strand`.
impl Display for Strand {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        match &self {
            Strand::Positive => write!(f, "+"),
            Strand::Negative => write!(f, "-"),
        }
    }
}
