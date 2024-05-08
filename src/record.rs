use crate::{Error, ErrorKind, Result};
use std::{
    fmt::{Display, Formatter},
    path::PathBuf,
    str::FromStr,
};

pub enum Record {
    Protein(ProteinRecord),
    Dna(DNARecord),
}

/// A record which could either be from a protein search or a
/// DNA search. As such, some methods will return `None` if the
/// record is from a protein search, and vice versa.
impl Record {
    /// The name of the target sequence or profile.
    pub fn target_name(&self) -> String {
        match self {
            Record::Protein(record) => record.target_name(),
            Record::Dna(record) => record.target_name(),
        }
    }
    /// The accession of the target sequence or profile, or ’-’ if none.
    pub fn target_accession(&self) -> String {
        match self {
            Record::Protein(record) => record.target_accession(),
            Record::Dna(record) => record.target_accession(),
        }
    }
    /// The name of the query sequence or profile.
    pub fn query_name(&self) -> String {
        match self {
            Record::Protein(record) => record.query_name(),
            Record::Dna(record) => record.query_name(),
        }
    }
    /// The accession of the query sequence or profile, or '-' if none.
    pub fn query_accession(&self) -> String {
        match self {
            Record::Protein(record) => record.query_accession(),
            Record::Dna(record) => record.query_accession(),
        }
    }
    /// The expectation value (statistical significance) of the target.
    /// This is a per query E-value; i.e. calculated as the expected
    /// number of false positives achieving this comparison’s score for
    /// a single query against the Z sequences in the target dataset.
    /// If you search with multiple queries and if you want to control
    /// the overall false positive rate of that search rather than the false
    /// positive rate per query, you will want to multiply this per query
    /// E-value by how many queries you’re doing. Protein (like) records only.
    pub fn e_value_full(&self) -> Option<f32> {
        match self {
            Record::Protein(record) => Some(record.e_value_full()),
            Record::Dna(_) => None,
        }
    }
    /// The score (in bits) for this target/query comparison. It includes
    /// the biased-composition correction (the “null2” model).
    pub fn score_full(&self) -> Option<f32> {
        match self {
            Record::Protein(record) => Some(record.score_full()),
            Record::Dna(record) => Some(record.score()),
        }
    }
    /// The biased-composition correction: the bit score difference
    /// contributed by the null2 model. High bias scores may be a
    /// red flag for a false positive, especially when the bias score is as
    /// large or larger than the overall bit score. It is difficult to
    /// correct for all possible ways in which a nonrandom but non-homologous
    /// biological sequences can appear to be similar, such as short-period
    /// tandem repeats, so there are cases where the bias correction is not
    /// strong enough (creating false positives). Protein (like) records only.
    pub fn bias_full(&self) -> Option<f32> {
        match self {
            Record::Protein(record) => Some(record.bias_full()),
            Record::Dna(record) => Some(record.bias()),
        }
    }
    /// The E-value if only the single best-scoring domain envelope were
    /// found in the sequence, and none of the others. If this E-value isn’t
    /// good, but the full sequence E-value is good, this is a potential
    /// red flag. Weak hits, none of which are good enough on their own,
    /// are summing up to lift the sequence up to a high score. Whether this is Good
    /// or Bad is not clear; the sequence may contain several weak homologous
    /// domains, or it might contain a repetitive sequence that is hitting by
    /// chance (i.e. once one repeat hits, all the repeats hit). Protein (like)
    /// records only.
    pub fn e_value_best(&self) -> Option<f32> {
        match self {
            Record::Protein(record) => Some(record.e_value_best()),
            Record::Dna(_) => None,
        }
    }
    /// The bit score if only the single best-scoring domain
    /// envelope were found in the sequence, and none of the others.
    /// (Inclusive of the null2 bias correction). Protein (like) records only.
    pub fn score_best(&self) -> Option<f32> {
        match self {
            Record::Protein(record) => Some(record.score_best()),
            Record::Dna(_) => None,
        }
    }
    /// The null2 bias correction that was applied to the bit score
    /// of the single best-scoring domain. Protein (like) records only.
    pub fn bias_best(&self) -> Option<f32> {
        match self {
            Record::Protein(record) => Some(record.bias_best()),
            Record::Dna(_) => None,
        }
    }
    /// Expected number of domains, as calculated by posterior decoding on
    /// the mean number of begin states used in the alignment ensemble.
    /// Protein (like) records only.
    pub fn exp(&self) -> Option<f32> {
        match self {
            Record::Protein(record) => Some(record.exp()),
            Record::Dna(_) => None,
        }
    }

    /// Number of discrete regions defined, as calculated by heuristics
    /// applied to posterior decoding of begin/end state positions
    /// in the alignment ensemble. The number of regions will generally
    /// be close to the expected number of domains. The more different
    /// the two numbers are, the less discrete the regions appear to be, in
    /// terms of probability mass. This usually means one of two things.
    /// On the one hand, weak homologous domains may be difficult for
    /// the heuristics to identify clearly. On the other hand, repetitive se-
    /// quence may appear to have a high expected domain number (from
    /// lots of crappy possible alignments in the ensemble, no one of
    /// which is very convincing on its own, so no one region is discretely
    /// well-defined). Protein (like) records only.
    pub fn reg(&self) -> Option<i32> {
        match self {
            Record::Protein(record) => Some(record.reg()),
            Record::Dna(_) => None,
        }
    }

    /// Number of regions that appeared to be multidomain, and therefore were
    /// passed to stochastic traceback clustering for further resolution down
    /// to one or more envelopes. This number is often zero. Protein (like) records only.
    pub fn clu(&self) -> Option<i32> {
        match self {
            Record::Protein(record) => Some(record.clu()),
            Record::Dna(_) => None,
        }
    }

    /// For envelopes that were defined by stochastic traceback
    /// clustering, how many of them overlap other envelopes.
    /// Protein (like) records only.
    pub fn ov(&self) -> Option<i32> {
        match self {
            Record::Protein(record) => Some(record.ov()),
            Record::Dna(_) => None,
        }
    }

    /// The total number of envelopes defined, both by single envelope regions
    /// and by stochastic traceback clustering into one or more envelopes per region.
    pub fn env(&self) -> Option<i32> {
        match self {
            Record::Protein(record) => Some(record.env()),
            Record::Dna(_) => None,
        }
    }

    /// Number of domains defined. In general, this is the same as the number of
    /// envelopes: for each envelope, we find an MEA (maximum expected accuracy)
    /// alignment, which defines the endpoints of the alignable domain. Protein
    /// (like) records only.
    pub fn dom(&self) -> Option<i32> {
        match self {
            Record::Protein(record) => Some(record.dom()),
            Record::Dna(_) => None,
        }
    }

    /// Number of domains satisfying reporting thresholds. If you’ve also saved a
    /// `--domtblout` file, there will be one line in it for each reported domain.
    /// Protein (like) records only.
    pub fn rep(&self) -> Option<i32> {
        match self {
            Record::Protein(record) => Some(record.rep()),
            Record::Dna(_) => None,
        }
    }

    /// Number of domains satisying inclusion thresholds. Protein (like)
    /// records only.
    pub fn inc(&self) -> Option<i32> {
        match self {
            Record::Protein(record) => Some(record.inc()),
            Record::Dna(_) => None,
        }
    }

    /// The position in the hmm at which the hit starts. DNA records only.
    pub fn hmm_from(&self) -> Option<i32> {
        match self {
            Record::Protein(_) => None,
            Record::Dna(record) => Some(record.hmm_from()),
        }
    }

    /// The position in the hmm at which the hit ends. DNA records only.
    pub fn hmm_to(&self) -> Option<i32> {
        match self {
            Record::Protein(_) => None,
            Record::Dna(record) => Some(record.hmm_to()),
        }
    }

    /// The position in the target sequence at which the hit starts.
    /// DNA records only.
    pub fn ali_from(&self) -> Option<i32> {
        match self {
            Record::Protein(_) => None,
            Record::Dna(record) => Some(record.ali_from()),
        }
    }
    /// The position in the target sequence at which the hit ends. DNA
    /// records only.
    pub fn ali_to(&self) -> Option<i32> {
        match self {
            Record::Protein(_) => None,
            Record::Dna(record) => Some(record.ali_to()),
        }
    }
    /// The position in the target sequence where the surrounding envelope starts.
    /// DNA records only.
    pub fn env_from(&self) -> Option<i32> {
        match self {
            Record::Protein(_) => None,
            Record::Dna(record) => Some(record.env_from()),
        }
    }
    /// The position in the target sequence at which the surrounding envelope ends.
    /// DNA records only.
    pub fn env_to(&self) -> Option<i32> {
        match self {
            Record::Protein(_) => None,
            Record::Dna(record) => Some(record.env_to()),
        }
    }
    /// The length of the target sequence. DNA records only.
    pub fn sq_len(&self) -> Option<i32> {
        match self {
            Record::Protein(_) => None,
            Record::Dna(record) => Some(record.sq_len()),
        }
    }
    /// The strand on which the hit was found
    /// (“-" when alifrom>ali to). DNA records only.
    pub fn strand(&self) -> Option<Strand> {
        match self {
            Record::Protein(_) => None,
            Record::Dna(record) => Some(record.strand()),
        }
    }
    /// The expectation value (statistical significance) of the target
    /// as above.
    pub fn e_value(&self) -> Option<f32> {
        match self {
            Record::Protein(_) => None,
            Record::Dna(record) => Some(record.e_value()),
        }
    }
}

#[derive(Default, PartialEq, Eq, Debug, Clone, Copy)]
pub enum Program {
    #[default]
    None,
    Nhmmer,    // test done
    Nhmmscan,  // test done
    Jackhmmer, // test done
    Hmmscan,   // test done
    Hmmsearch,
    Phmmer, // test done
}

impl FromStr for Program {
    type Err = Error;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s {
            "nhmmer" => Ok(Program::Nhmmer),
            "nhmmscan" => Ok(Program::Nhmmscan),
            "jackhmmer" => Ok(Program::Jackhmmer),
            "hmmscan" => Ok(Program::Hmmscan),
            "hmmsearch" => Ok(Program::Hmmsearch),
            "phmmer" => Ok(Program::Phmmer),
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
        self.program
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

#[derive(Debug)]
pub struct ProteinRecord {
    target_name: String,
    target_accession: String,
    query_name: String,
    query_accession: String,
    e_value_full: f32,
    score_full: f32,
    bias_full: f32,
    e_value_best: f32,
    score_best: f32,
    bias_best: f32,
    exp: f32,
    reg: i32,
    clu: i32,
    ov: i32,
    env: i32,
    dom: i32,
    rep: i32,
    inc: i32,
}

impl ProteinRecord {
    pub fn new(
        target_name: String,
        target_accession: String,
        query_name: String,
        query_accession: String,
        e_value_full: f32,
        score_full: f32,
        bias_full: f32,
        e_value_best: f32,
        score_best: f32,
        bias_best: f32,
        exp: f32,
        reg: i32,
        clu: i32,
        ov: i32,
        env: i32,
        dom: i32,
        rep: i32,
        inc: i32,
    ) -> Self {
        ProteinRecord {
            target_name,
            target_accession,
            query_name,
            query_accession,
            e_value_full,
            score_full,
            bias_full,
            e_value_best,
            score_best,
            bias_best,
            exp,
            reg,
            clu,
            ov,
            env,
            dom,
            rep,
            inc,
        }
    }

    pub fn target_name(&self) -> String {
        self.target_name.clone()
    }

    pub fn target_accession(&self) -> String {
        self.target_accession.clone()
    }

    pub fn query_name(&self) -> String {
        self.query_name.clone()
    }

    pub fn query_accession(&self) -> String {
        self.query_accession.clone()
    }

    pub fn e_value_full(&self) -> f32 {
        self.e_value_full
    }

    pub fn score_full(&self) -> f32 {
        self.score_full
    }

    pub fn bias_full(&self) -> f32 {
        self.bias_full
    }

    pub fn e_value_best(&self) -> f32 {
        self.e_value_best
    }

    pub fn score_best(&self) -> f32 {
        self.score_best
    }

    pub fn bias_best(&self) -> f32 {
        self.bias_best
    }

    pub fn exp(&self) -> f32 {
        self.exp
    }

    pub fn reg(&self) -> i32 {
        self.reg
    }

    pub fn clu(&self) -> i32 {
        self.clu
    }

    pub fn ov(&self) -> i32 {
        self.ov
    }

    pub fn env(&self) -> i32 {
        self.env
    }

    pub fn dom(&self) -> i32 {
        self.dom
    }

    pub fn rep(&self) -> i32 {
        self.rep
    }

    pub fn inc(&self) -> i32 {
        self.inc
    }

    pub fn set_target_name(&mut self, target_name: String) {
        self.target_name = target_name;
    }

    pub fn set_target_accession(&mut self, target_accession: String) {
        self.target_accession = target_accession;
    }

    pub fn set_query_name(&mut self, query_name: String) {
        self.query_name = query_name;
    }

    pub fn set_query_accession(&mut self, query_accession: String) {
        self.query_accession = query_accession;
    }

    pub fn set_e_value_full(&mut self, e_value_full: f32) {
        self.e_value_full = e_value_full;
    }

    pub fn set_score_full(&mut self, score_full: f32) {
        self.score_full = score_full;
    }

    pub fn set_bias_full(&mut self, bias_full: f32) {
        self.bias_full = bias_full;
    }

    pub fn set_e_value_best(&mut self, e_value_best: f32) {
        self.e_value_best = e_value_best;
    }

    pub fn set_score_best(&mut self, score_best: f32) {
        self.score_best = score_best;
    }

    pub fn set_bias_best(&mut self, bias_best: f32) {
        self.bias_best = bias_best;
    }

    pub fn set_exp(&mut self, exp: f32) {
        self.exp = exp;
    }

    pub fn set_reg(&mut self, reg: i32) {
        self.reg = reg;
    }

    pub fn set_clu(&mut self, clu: i32) {
        self.clu = clu;
    }

    pub fn set_ov(&mut self, ov: i32) {
        self.ov = ov;
    }

    pub fn set_env(&mut self, env: i32) {
        self.env = env;
    }

    pub fn set_dom(&mut self, dom: i32) {
        self.dom = dom;
    }

    pub fn set_rep(&mut self, rep: i32) {
        self.rep = rep;
    }

    pub fn set_inc(&mut self, inc: i32) {
        self.inc = inc;
    }
}

/// A record in a HMMER tblout file
/// specific to DNA related searches.
#[derive(Debug)]
pub struct DNARecord {
    target_name: String,
    target_accession: String,
    query_name: String,
    query_accession: String,
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

impl DNARecord {
    pub fn new(
        target_name: String,
        target_accession: String,
        query_name: String,
        query_accession: String,
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
        DNARecord {
            target_name,
            target_accession,
            query_name,
            query_accession,
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

    pub fn target_accession(&self) -> String {
        self.target_accession.clone()
    }

    pub fn query_name(&self) -> String {
        self.query_name.clone()
    }

    pub fn query_accession(&self) -> String {
        self.query_accession.clone()
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
