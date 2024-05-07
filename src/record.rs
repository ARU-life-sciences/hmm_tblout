use crate::{Error, ErrorKind, Result};
use std::{
    fmt::{Display, Formatter},
    path::PathBuf,
    str::FromStr,
};

#[derive(Default, PartialEq, Eq, Debug, Clone, Copy)]
pub enum Program {
    #[default]
    None,
    Nhmmer,
}

impl FromStr for Program {
    type Err = Error;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s {
            "nhmmer" => Ok(Program::Nhmmer),
            _ => Err(Error::new(ErrorKind::Meta(format!(
                "The program \"{}\" is not supported.",
                s
            )))),
        }
    }
}

#[derive(Default)]
pub struct Meta {
    program: Program,
    version: String,
    pipeline_mode: String,
    query_file: PathBuf,
    target_file: PathBuf,
    options: String,
    current_dir: PathBuf,
    date: String,
}

impl Meta {
    /// Get the program information
    pub fn program(&self) -> Program {
        self.program.clone()
    }

    /// Set program
    pub fn set_program(&mut self, program: Program) {
        self.program = program;
    }

    /// Get the version
    pub fn version(&self) -> String {
        self.version.clone()
    }

    /// Set version
    pub fn set_version(&mut self, version: String) {
        self.version = version;
    }

    /// Get the pipeline mode
    pub fn pipeline_mode(&self) -> String {
        self.pipeline_mode.clone()
    }

    /// Set pipeline mode
    pub fn set_pipeline_mode(&mut self, pipeline_mode: String) {
        self.pipeline_mode = pipeline_mode;
    }

    /// Get the path to the query file
    pub fn query_file(&self) -> PathBuf {
        self.query_file.clone()
    }

    /// Set path
    pub fn set_query_file(&mut self, query_file: PathBuf) {
        self.query_file = query_file;
    }

    /// Get the path to the target file
    pub fn target_file(&self) -> PathBuf {
        self.target_file.clone()
    }

    /// Set target file
    pub fn set_target_file(&mut self, target_file: PathBuf) {
        self.target_file = target_file;
    }

    /// Get the options used to run the program
    pub fn options(&self) -> String {
        self.options.clone()
    }

    /// Set options
    pub fn set_options(&mut self, options: String) {
        self.options = options;
    }

    /// Get the current directory
    pub fn current_dir(&self) -> PathBuf {
        self.current_dir.clone()
    }

    /// Set current directory
    pub fn set_current_dir(&mut self, current_dir: PathBuf) {
        self.current_dir = current_dir;
    }

    /// Get the date the program was run
    pub fn date(&self) -> String {
        self.date.clone()
    }

    /// Set the date
    pub fn set_date(&mut self, date: String) {
        self.date = date;
    }
}

/// A record in a HMMER tblout file
#[derive(Debug)]
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
