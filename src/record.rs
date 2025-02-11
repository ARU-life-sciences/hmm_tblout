use crate::{Error, ErrorKind, Result};
use std::{
    fmt::{Display, Formatter},
    path::PathBuf,
    str::FromStr,
};

/// A record in a HMMER tblout file. Can either be a protein
/// record or a DNA record.
pub enum Record {
    /// A protein record.
    Protein(ProteinRecord),
    /// A DNA record.
    Dna(DNARecord),
    /// A cmsearch/cmscan record.
    CM(CMRecord),
}

/// A record which could either be from a protein search or a
/// DNA search. As such, some methods will return `None` if the
/// record is from a protein search, and vice versa. The descriptions
/// from the HMMER user guide are used here (v3.4 August 2023).
impl Record {
    /// The name of the target sequence or profile.
    pub fn target_name(&self) -> String {
        match self {
            Record::Protein(record) => record.target_name(),
            Record::Dna(record) => record.target_name(),
            Record::CM(cmrecord) => cmrecord.target_name(),
        }
    }
    /// The accession of the target sequence or profile, or ’-’ if none.
    pub fn target_accession(&self) -> String {
        match self {
            Record::Protein(record) => record.target_accession(),
            Record::Dna(record) => record.target_accession(),
            Record::CM(cmrecord) => cmrecord.target_accession(),
        }
    }
    /// The name of the query sequence or profile.
    pub fn query_name(&self) -> String {
        match self {
            Record::Protein(record) => record.query_name(),
            Record::Dna(record) => record.query_name(),
            Record::CM(cmrecord) => cmrecord.query_name(),
        }
    }
    /// The accession of the query sequence or profile, or '-' if none.
    pub fn query_accession(&self) -> String {
        match self {
            Record::Protein(record) => record.query_accession(),
            Record::Dna(record) => record.query_accession(),
            Record::CM(cmrecord) => cmrecord.query_accession(),
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
            Record::CM(_) => None,
        }
    }
    /// The score (in bits) for this target/query comparison. It includes
    /// the biased-composition correction (the “null2” model). Protein (like)
    /// records only.
    pub fn score_full(&self) -> Option<f32> {
        match self {
            Record::Protein(record) => Some(record.score_full()),
            Record::Dna(_) => None,
            Record::CM(_) => None,
        }
    }

    /// The score (in bits) for this target/query comparison. It includes the biased-composition correction (the
    /// “null3” model for CM hits, or the “null2” model for HMM hits). DNA and cmsearch/cmscan.
    pub fn score(&self) -> Option<f32> {
        match self {
            Record::Protein(_) => None,
            Record::Dna(record) => Some(record.score()),
            Record::CM(cmrecord) => Some(cmrecord.score()),
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
            Record::Dna(_) => None,
            Record::CM(_) => None,
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
            Record::CM(_) => None,
        }
    }
    /// The bit score if only the single best-scoring domain
    /// envelope were found in the sequence, and none of the others.
    /// (Inclusive of the null2 bias correction). Protein (like) records only.
    pub fn score_best(&self) -> Option<f32> {
        match self {
            Record::Protein(record) => Some(record.score_best()),
            Record::Dna(_) => None,
            Record::CM(_) => None,
        }
    }
    /// The null2 bias correction that was applied to the bit score
    /// of the single best-scoring domain. Protein (like) records only.
    pub fn bias_best(&self) -> Option<f32> {
        match self {
            Record::Protein(record) => Some(record.bias_best()),
            Record::Dna(_) => None,
            Record::CM(_) => None,
        }
    }

    /// The biased-composition correction: the bit score difference contributed by the null3 model for CM hits, or
    /// the null2 model for HMM hits. High bias scores may be a red flag for a false positive. It is difficult to correct for all
    /// possible ways in which a nonrandom but nonhomologous biological sequences can appear to be similar, such as
    /// short-period tandem repeats, so there are cases where the bias correction is not strong enough (creating false
    /// positives). DNA and cmsearch/cmscan only.
    pub fn bias(&self) -> Option<f32> {
        match self {
            Record::Protein(_) => None,
            Record::Dna(record) => Some(record.bias()),
            Record::CM(cmrecord) => Some(cmrecord.bias()),
        }
    }

    /// Expected number of domains, as calculated by posterior decoding on
    /// the mean number of begin states used in the alignment ensemble.
    /// Protein (like) records only.
    pub fn exp(&self) -> Option<f32> {
        match self {
            Record::Protein(record) => Some(record.exp()),
            Record::Dna(_) => None,
            Record::CM(_) => None,
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
            Record::CM(_) => None,
        }
    }

    /// Number of regions that appeared to be multidomain, and therefore were
    /// passed to stochastic traceback clustering for further resolution down
    /// to one or more envelopes. This number is often zero. Protein (like) records only.
    pub fn clu(&self) -> Option<i32> {
        match self {
            Record::Protein(record) => Some(record.clu()),
            Record::Dna(_) => None,
            Record::CM(_) => None,
        }
    }

    /// For envelopes that were defined by stochastic traceback
    /// clustering, how many of them overlap other envelopes.
    /// Protein (like) records only.
    pub fn ov(&self) -> Option<i32> {
        match self {
            Record::Protein(record) => Some(record.ov()),
            Record::Dna(_) => None,
            Record::CM(_) => None,
        }
    }

    /// The total number of envelopes defined, both by single envelope regions
    /// and by stochastic traceback clustering into one or more envelopes per region.
    /// Protein (like) records only.
    pub fn env(&self) -> Option<i32> {
        match self {
            Record::Protein(record) => Some(record.env()),
            Record::Dna(_) => None,
            Record::CM(_) => None,
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
            Record::CM(_) => None,
        }
    }

    /// Number of domains satisfying reporting thresholds. If you’ve also saved a
    /// `--domtblout` file, there will be one line in it for each reported domain.
    /// Protein (like) records only.
    pub fn rep(&self) -> Option<i32> {
        match self {
            Record::Protein(record) => Some(record.rep()),
            Record::Dna(_) => None,
            Record::CM(_) => None,
        }
    }

    /// Number of domains satisying inclusion thresholds. Protein (like)
    /// records only.
    pub fn inc(&self) -> Option<i32> {
        match self {
            Record::Protein(record) => Some(record.inc()),
            Record::Dna(_) => None,
            Record::CM(_) => None,
        }
    }

    /// Indicates whether or not this hit achieves the inclusion threshold: ’!’ if it does, ’?’ if it does not (and rather
    /// only achieves the reporting threshold). By default, the inclusion threshold is an E-value of 0.01 and the reporting
    /// threshold is an E-value of 10.0, but these can be changed with command line options as described in the manual
    /// pages. cmsearch/cmscan only.
    pub fn cm_inc(&self) -> Option<char> {
        match self {
            Record::Protein(_) => None,
            Record::Dna(_) => None,
            Record::CM(cmrecord) => Some(cmrecord.inc()),
        }
    }

    /// The position in the hmm at which the hit starts. DNA records only.
    pub fn hmm_from(&self) -> Option<i32> {
        match self {
            Record::Protein(_) => None,
            Record::Dna(record) => Some(record.hmm_from()),
            Record::CM(_) => None,
        }
    }

    /// The position in the hmm at which the hit ends. DNA records only.
    pub fn hmm_to(&self) -> Option<i32> {
        match self {
            Record::Protein(_) => None,
            Record::Dna(record) => Some(record.hmm_to()),
            Record::CM(_) => None,
        }
    }

    /// The position in the target sequence at which the hit starts.
    /// DNA records only.
    pub fn ali_from(&self) -> Option<i32> {
        match self {
            Record::Protein(_) => None,
            Record::Dna(record) => Some(record.ali_from()),
            Record::CM(_) => None,
        }
    }

    /// The position in the target sequence at which the hit ends. DNA
    /// records only.
    pub fn ali_to(&self) -> Option<i32> {
        match self {
            Record::Protein(_) => None,
            Record::Dna(record) => Some(record.ali_to()),
            Record::CM(_) => None,
        }
    }

    /// The position in the target sequence where the surrounding envelope starts.
    /// DNA records only.
    pub fn env_from(&self) -> Option<i32> {
        match self {
            Record::Protein(_) => None,
            Record::Dna(record) => Some(record.env_from()),
            Record::CM(_) => None,
        }
    }

    /// The position in the target sequence at which the surrounding envelope ends.
    /// DNA records only.
    pub fn env_to(&self) -> Option<i32> {
        match self {
            Record::Protein(_) => None,
            Record::Dna(record) => Some(record.env_to()),
            Record::CM(_) => None,
        }
    }

    /// The length of the target sequence. DNA records only.
    pub fn sq_len(&self) -> Option<i32> {
        match self {
            Record::Protein(_) => None,
            Record::Dna(record) => Some(record.sq_len()),
            Record::CM(_) => None,
        }
    }

    /// The strand on which the hit was found
    /// (“-" when alifrom>ali to). DNA and cmsearch/cmscan records only.
    pub fn strand(&self) -> Option<Strand> {
        match self {
            Record::Protein(_) => None,
            Record::Dna(record) => Some(record.strand()),
            Record::CM(cmrecord) => Some(cmrecord.strand()),
        }
    }

    /// The expectation value (statistical significance) of the target
    /// as above. DNA and cmsearch/cmscan records only.
    pub fn e_value(&self) -> Option<f32> {
        match self {
            Record::Protein(_) => None,
            Record::Dna(record) => Some(record.e_value()),
            Record::CM(cmrecord) => Some(cmrecord.e_value()),
        }
    }

    /// Which type of model was used to compute the final score. Either ’cm’ or ’hmm’. A CM is used
    // to compute the final hit scores unless the model has zero basepairs or the --hmmonly option is used, in which
    // case a HMM will be used. cmsearch/cmscan only.
    pub fn mdl(&self) -> Option<String> {
        match self {
            Record::Protein(_) => None,
            Record::Dna(_) => None,
            Record::CM(cmrecord) => Some(cmrecord.mdl()),
        }
    }

    /// The start of the alignment of this hit with respect to the profile (CM or HMM),
    /// numbered 1..N for a profile of N consensus positions. cmsearch/cmscan only.
    pub fn mdl_from(&self) -> Option<i32> {
        match self {
            Record::Protein(_) => None,
            Record::Dna(_) => None,
            Record::CM(cmrecord) => Some(cmrecord.mdl_from()),
        }
    }

    /// The end of the alignment of this hit with respect to the profile (CM or HMM),
    /// numbered 1..N for a profile of N consensus positions. cmsearch/cmscan only.
    pub fn mdl_to(&self) -> Option<i32> {
        match self {
            Record::Protein(_) => None,
            Record::Dna(_) => None,
            Record::CM(cmrecord) => Some(cmrecord.mdl_to()),
        }
    }

    /// The start of the alignment of this hit with respect to the sequence, numbered 1..L
    /// for a sequence of L residues. cmsearch/cmscan only.
    pub fn seq_from(&self) -> Option<i32> {
        match self {
            Record::Protein(_) => None,
            Record::Dna(_) => None,
            Record::CM(cmrecord) => Some(cmrecord.seq_from()),
        }
    }

    /// The end of the alignment of this hit with respect to the sequence, numbered 1..L for a
    /// sequence of L residues. cmsearch/cmscan only.
    pub fn seq_to(&self) -> Option<i32> {
        match self {
            Record::Protein(_) => None,
            Record::Dna(_) => None,
            Record::CM(cmrecord) => Some(cmrecord.seq_to()),
        }
    }

    /// Indicates if this is predicted to be a truncated CM hit or not. This will be “no” if it is a CM hit that is not
    /// predicted to be truncated by the end of the sequence, “5’ ” or “3’ ” if the hit is predicted to have one or more 5’ or
    /// 3’ residues missing due to a artificial truncation of the sequence, or “5’&3”’ if the hit is predicted to have one or
    /// more 5’ residues missing and one or more 3’ residues missing. If the hit is an HMM hit, this will always be ’-’. cmsearch/cmscan only.
    pub fn trunc(&self) -> Option<String> {
        match self {
            Record::Protein(_) => None,
            Record::Dna(_) => None,
            Record::CM(cmrecord) => Some(cmrecord.trunc()),
        }
    }

    /// Indicates what “pass” of the pipeline the hit was detected on. This is probably only useful for testing and
    /// debugging. Non-truncated hits are found on the first pass, truncated hits are found on successive passes. cmsearch/cmscan only.
    pub fn pass(&self) -> Option<i32> {
        match self {
            Record::Protein(_) => None,
            Record::Dna(_) => None,
            Record::CM(cmrecord) => Some(cmrecord.pass()),
        }
    }

    /// Fraction of G and C nucleotides in the hit. cmsearch/cmscan only.
    pub fn gc(&self) -> Option<f32> {
        match self {
            Record::Protein(_) => None,
            Record::Dna(_) => None,
            Record::CM(cmrecord) => Some(cmrecord.gc()),
        }
    }

    /// A description, as free text. Available in DNA, protein, and cmsearch/cmscan records.
    pub fn description(&self) -> String {
        match self {
            Record::Protein(record) => record.description(),
            Record::Dna(record) => record.description(),
            Record::CM(cmrecord) => cmrecord.description(),
        }
    }
}

#[derive(Default, PartialEq, Eq, Debug, Clone, Copy)]
/// The program used to generate the output.
pub enum Program {
    #[default]
    /// The program is unknown. This is an error.
    None,
    /// The program used was `nhmmer`.
    Nhmmer, // test done
    /// The program used was `nhmmscan`.
    Nhmmscan, // test done
    /// The program used was `jackhmmer`.
    Jackhmmer, // test done
    /// The program used was `hmmscan`.
    Hmmscan, // test done
    /// The program used was `hmmsearch`.
    Hmmsearch,
    /// The program used was `phmmer`.
    Phmmer, // test done
    /// The program used was `cmsearch`.
    Cmsearch,
    /// The program used was `cmscan`.
    Cmscan,
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
            "cmsearch" => Ok(Program::Cmsearch),
            "cmscan" => Ok(Program::Cmscan),
            _ => Err(Error::new(ErrorKind::Meta(format!(
                "The program \"{}\" is not supported.",
                s
            )))),
        }
    }
}

impl Display for Program {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        let inner = match self {
            Program::Nhmmer => "nhmmer",
            Program::Nhmmscan => "nhmmscan",
            Program::Jackhmmer => "jackhmmer",
            Program::Hmmscan => "hmmscan",
            Program::Hmmsearch => "hmmsearch",
            Program::Phmmer => "phmmer",
            Program::Cmsearch => "cmsearch",
            Program::Cmscan => "cmscan",
            // FIXME: make this error better.
            Program::None => return Err(std::fmt::Error),
        };
        write!(f, "{}", inner)
    }
}

#[derive(Default, Clone)]
pub struct Header {
    protein_only: Option<String>,
    columns: String,
    dashes: String,
}

impl Display for Header {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        fn trim_newline(s: &mut String) {
            if s.ends_with('\n') {
                s.pop();
                if s.ends_with('\r') {
                    s.pop();
                }
            }
        }
        // if the protein_only field is Some, print it.
        if let Some(mut p) = self.protein_only.clone() {
            trim_newline(&mut p);
            writeln!(f, "{}", p)?;
        }

        let mut cols = self.columns.clone();
        let mut dashes = self.dashes.clone();
        trim_newline(&mut cols);
        trim_newline(&mut dashes);

        writeln!(f, "{}", cols)?;
        // should this also be a newline?
        write!(f, "{}", dashes)
    }
}

impl Header {
    pub fn new(protein_only: Option<String>, columns: String, dashes: String) -> Self {
        Header {
            protein_only,
            columns,
            dashes,
        }
    }
    // calculate a vector of column widths
    pub fn calculate_dashes(&self) -> Vec<usize> {
        let mut lengths = Vec::new();
        let mut current_length = 0;
        let mut space_count = 0;
        let mut in_dash_run = false;

        for c in self.dashes.chars() {
            match c {
                '-' | '#' => {
                    if space_count > 1 {
                        current_length += space_count - 1; // Include (n-1) spaces
                    }
                    space_count = 0;
                    current_length += 1;
                    in_dash_run = true;
                }
                ' ' => {
                    if in_dash_run {
                        if space_count == 0 {
                            lengths.push(current_length); // Push previous dash run
                            current_length = 0;
                        }
                        space_count += 1;
                    }
                }
                _ => {
                    if current_length > 0 {
                        lengths.push(current_length);
                        current_length = 0;
                    }
                    space_count = 0;
                    in_dash_run = false;
                }
            }
        }

        if current_length > 0 {
            lengths.push(current_length);
        }

        lengths
    }

    /// Get the protein_only field.
    pub fn get_protein_only(&self) -> Option<String> {
        self.protein_only.clone()
    }

    /// Set the protein_only field.
    pub fn set_protein_only(&mut self, protein_only: String) {
        self.protein_only = Some(protein_only);
    }

    /// Get the columns field.
    pub fn get_columns(&self) -> String {
        self.columns.clone()
    }

    /// Set the columns field.
    pub fn set_columns(&mut self, columns: String) {
        self.columns = columns;
    }

    /// Get the dashes field.
    pub fn get_dashes(&self) -> String {
        self.dashes.clone()
    }

    /// Set the dashes field.
    pub fn set_dashes(&mut self, dashes: String) {
        self.dashes = dashes;
    }
}

#[derive(Default, Clone)]
/// Metadata about the search that produced the HMMER tblout file.
pub struct Meta {
    /// The program used to generate the output.
    program: Program,
    /// The version of the program used to generate the output.
    version: String,
    /// The pipeline mode used to generate the output.
    pipeline_mode: String,
    /// The path to the query file.
    query_file: PathBuf,
    /// The path to the target file.
    target_file: PathBuf,
    /// The options used to run the program.
    options: String,
    /// The current directory.
    current_dir: PathBuf,
    /// The date the program was run.
    date: String,
}

// Implement Display for Meta.
impl Display for Meta {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        writeln!(f, "#")?;
        writeln!(f, "# Program:         {}", self.program)?;
        writeln!(f, "# Version:         {}", self.version)?;
        writeln!(f, "# Pipeline mode:   {}", self.pipeline_mode)?;
        writeln!(f, "# Query file:      {}", self.query_file.display())?;
        writeln!(f, "# Target file:     {}", self.target_file.display())?;
        writeln!(f, "# Options:         {}", self.options)?;
        writeln!(f, "# Current dir:     {}", self.current_dir.display())?;
        writeln!(f, "# Date:            {}", self.date)?;
        write!(f, "# [ok]")
    }
}

impl Meta {
    /// Get the program information.
    pub fn program(&self) -> Program {
        self.program
    }

    /// Set program.
    pub fn set_program(&mut self, program: Program) {
        self.program = program;
    }

    /// Get the version.
    pub fn version(&self) -> String {
        self.version.clone()
    }

    /// Set version.
    pub fn set_version(&mut self, version: String) {
        self.version = version;
    }

    /// Get the pipeline mode.
    pub fn pipeline_mode(&self) -> String {
        self.pipeline_mode.clone()
    }

    /// Set pipeline mode.
    pub fn set_pipeline_mode(&mut self, pipeline_mode: String) {
        self.pipeline_mode = pipeline_mode;
    }

    /// Get the path to the query file.
    pub fn query_file(&self) -> PathBuf {
        self.query_file.clone()
    }

    /// Set path.
    pub fn set_query_file(&mut self, query_file: PathBuf) {
        self.query_file = query_file;
    }

    /// Get the path to the target file.
    pub fn target_file(&self) -> PathBuf {
        self.target_file.clone()
    }

    /// Set target file.
    pub fn set_target_file(&mut self, target_file: PathBuf) {
        self.target_file = target_file;
    }

    /// Get the options used to run the program.
    pub fn options(&self) -> String {
        self.options.clone()
    }

    /// Set options.
    pub fn set_options(&mut self, options: String) {
        self.options = options;
    }

    /// Get the current directory.
    pub fn current_dir(&self) -> PathBuf {
        self.current_dir.clone()
    }

    /// Set current directory.
    pub fn set_current_dir(&mut self, current_dir: PathBuf) {
        self.current_dir = current_dir;
    }

    /// Get the date the program was run.
    pub fn date(&self) -> String {
        self.date.clone()
    }

    /// Set the date.
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
    description: String,
    col_sizes: Vec<usize>,
}

impl Display for ProteinRecord {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        write!(
            f,
            "{:<width0$} {:<width1$} {:<width2$} {:<width3$} {:>width4$e} {:>width5$} {:>width6$} {:>width7$e} {:>width8$} {:>width9$} {:>width10$.1} {:>width11$} {:>width12$} {:>width13$} {:>width14$} {:>width15$} {:>width16$} {:>width17$} {:>width18$}",
            self.target_name,
            self.target_accession,
            self.query_name,
            self.query_accession,
            self.e_value_full,
            self.score_full,
            self.bias_full,
            self.e_value_best,
            self.score_best,
            self.bias_best,
            self.exp,
            self.reg,
            self.clu,
            self.ov,
            self.env,
            self.dom,
            self.rep,
            self.inc,
            self.description,
            width0 = self.col_sizes[0],
            width1 = self.col_sizes[1],
            width2 = self.col_sizes[2],
            width3 = self.col_sizes[3],
            width4 = self.col_sizes[4],
            width5 = self.col_sizes[5],
            width6 = self.col_sizes[6],
            width7 = self.col_sizes[7],
            width8 = self.col_sizes[8],
            width9 = self.col_sizes[9],
            width10 = self.col_sizes[10],
            width11 = self.col_sizes[11],
            width12 = self.col_sizes[12],
            width13 = self.col_sizes[13],
            width14 = self.col_sizes[14],
            width15 = self.col_sizes[15],
            width16 = self.col_sizes[16],
            width17 = self.col_sizes[17],
            width18 = self.col_sizes[18]
        )
    }
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
        description: String,
        col_sizes: Vec<usize>,
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
            description,
            col_sizes,
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

    pub fn description(&self) -> String {
        self.description.clone()
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

    pub fn set_description(&mut self, description: String) {
        self.description = description;
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
    description: String,
    col_sizes: Vec<usize>,
}

// a display implementation for `DNARecord`
// the original HMMER output is:
// "...deliberately space-
// delimited (rather than tab-delimited)
// and justified into aligned columns,
// so these files are suitable both for
// automated parsing and for human
// examination. I feel that tab-delimited
// data files are difficult for humans
// to examine and spot check. For this
// reason, I think tab-delimited files are
// a minor evil in the world. Although I
// occasionally receive shrieks of outrage
// about this, I still stubbornly feel that
// space-delimited files are just as easily
// parsed as tab-delimited files"

impl Display for DNARecord {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        write!(
            f,
            // tabs would be just as good, Eddy.
            "{:<width0$} {:<width1$} {:<width2$} {:<width3$} {:>width4$} {:>width5$} {:>width6$} {:>width7$} {:>width8$} {:>width9$} {:>width10$} {:^width11$} {:>width12$e} {:>width13$} {:>width14$} {:<width15$}",
            self.target_name,
            self.target_accession,
            self.query_name,
            self.query_accession,
            self.hmm_from,
            self.hmm_to,
            self.ali_from,
            self.ali_to,
            self.env_from,
            self.env_to,
            self.sq_len,
            self.strand.to_string(),
            self.e_value,
            self.score,
            self.bias,
            self.description,
            width0 = self.col_sizes[0],
            width1 = self.col_sizes[1],
            width2 = self.col_sizes[2],
            width3 = self.col_sizes[3],
            width4 = self.col_sizes[4],
            width5 = self.col_sizes[5],
            width6 = self.col_sizes[6],
            width7 = self.col_sizes[7],
            width8 = self.col_sizes[8],
            width9 = self.col_sizes[9],
            width10 = self.col_sizes[10],
            width11 = self.col_sizes[11],
            width12 = self.col_sizes[12],
            width13 = self.col_sizes[13],
            width14 = self.col_sizes[14],
            width15 = self.col_sizes[15],
        )
    }
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
        description: String,
        col_sizes: Vec<usize>,
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
            description,
            col_sizes,
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

    pub fn description(&self) -> String {
        self.description.clone()
    }
}

#[derive(Debug)]
/// The target hits tables produced from `cmsearch`
/// and `cmscan`. We cater *only* for format 1.
/// i.e. the default format.
pub struct CMRecord {
    target_name: String,
    target_accession: String,
    query_name: String,
    query_accession: String,
    mdl: String,
    mdl_from: i32,
    mdl_to: i32,
    seq_from: i32,
    seq_to: i32,
    strand: Strand,
    trunc: String,
    pass: i32,
    gc: f32,
    bias: f32,
    score: f32,
    e_value: f32,
    inc: char,
    description: String,
    col_sizes: Vec<usize>,
}

impl Display for CMRecord {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        write!(
            f,
            "{:<width0$} {:<width1$} {:<width2$} {:<width3$} {:<width4$} {:>width5$} {:>width6$} {:>width7$} {:>width8$} {:>width9$} {:>width10$} {:>width11$} {:>width12$} {:>width13$} {:>width14$} {:>width15$} {:^width16$} {:<width17$}",
            self.target_name,
            self.target_accession,
            self.query_name,
            self.query_accession,
            self.mdl,
            self.mdl_from,
            self.mdl_to,
            self.seq_from,
            self.seq_to,
            self.strand.to_string(),
            self.trunc,
            self.pass,
            self.gc,
            self.bias,
            self.score,
            self.e_value,
            self.inc,
            self.description,
            width0 = self.col_sizes[0],
            width1 = self.col_sizes[1],
            width2 = self.col_sizes[2],
            width3 = self.col_sizes[3],
            width4 = self.col_sizes[4],
            width5 = self.col_sizes[5],
            width6 = self.col_sizes[6],
            width7 = self.col_sizes[7],
            width8 = self.col_sizes[8],
            width9 = self.col_sizes[9],
            width10 = self.col_sizes[10],
            width11 = self.col_sizes[11],
            width12 = self.col_sizes[12],
            width13 = self.col_sizes[13],
            width14 = self.col_sizes[14],
            width15 = self.col_sizes[15],
            width16 = self.col_sizes[16],
            width17 = self.col_sizes[17]
        )
    }
}

impl CMRecord {
    pub fn new(
        target_name: String,
        target_accession: String,
        query_name: String,
        query_accession: String,
        mdl: String,
        mdl_from: i32,
        mdl_to: i32,
        seq_from: i32,
        seq_to: i32,
        strand: Strand,
        trunc: String,
        pass: i32,
        gc: f32,
        bias: f32,
        score: f32,
        e_value: f32,
        inc: char,
        description: String,
        col_sizes: Vec<usize>,
    ) -> Self {
        CMRecord {
            target_name,
            target_accession,
            query_name,
            query_accession,
            mdl,
            mdl_from,
            mdl_to,
            seq_from,
            seq_to,
            strand,
            trunc,
            pass,
            gc,
            bias,
            score,
            e_value,
            inc,
            description,
            col_sizes,
        }
    }

    /// The name of the target sequence or profile.
    pub fn target_name(&self) -> String {
        self.target_name.clone()
    }

    /// The accession of the target sequence or profile, or ’-’ if none.
    pub fn target_accession(&self) -> String {
        self.target_accession.clone()
    }
    /// The name of the query sequence or profile.
    pub fn query_name(&self) -> String {
        self.query_name.clone()
    }
    /// The accession of the query sequence or profile, or ’-’ if none.
    pub fn query_accession(&self) -> String {
        self.query_accession.clone()
    }
    /// Which type of model was used to compute the final score. Either ’cm’ or ’hmm’. A CM is used
    // to compute the final hit scores unless the model has zero basepairs or the --hmmonly option is used, in which
    // case a HMM will be used.
    pub fn mdl(&self) -> String {
        self.mdl.clone()
    }
    /// The start of the alignment of this hit with respect to the profile (CM or HMM),
    /// numbered 1..N for a profile of N consensus positions.
    pub fn mdl_from(&self) -> i32 {
        self.mdl_from
    }
    /// The end of the alignment of this hit with respect to the profile (CM or HMM),
    /// numbered 1..N for a profile of N consensus positions.
    pub fn mdl_to(&self) -> i32 {
        self.mdl_to
    }
    /// The start of the alignment of this hit with respect to the sequence, numbered 1..L
    /// for a sequence of L residues.
    pub fn seq_from(&self) -> i32 {
        self.seq_from
    }
    /// The end of the alignment of this hit with respect to the sequence, numbered 1..L for a
    /// sequence of L residues.
    pub fn seq_to(&self) -> i32 {
        self.seq_to
    }
    /// The strand on which the hit occurs on the sequence. ’+’ if the hit is on the top (Watson) strand, ’-’ if
    /// the hit is on the bottom (Crick) strand. If on the top strand, the “seq from” value will be less than or equal to the
    /// “seq to” value, else it will be greater than or equal to it.
    pub fn strand(&self) -> Strand {
        self.strand
    }
    /// Indicates if this is predicted to be a truncated CM hit or not. This will be “no” if it is a CM hit that is not
    /// predicted to be truncated by the end of the sequence, “5’ ” or “3’ ” if the hit is predicted to have one or more 5’ or
    /// 3’ residues missing due to a artificial truncation of the sequence, or “5’&3”’ if the hit is predicted to have one or
    /// more 5’ residues missing and one or more 3’ residues missing. If the hit is an HMM hit, this will always be ’-’.
    pub fn trunc(&self) -> String {
        self.trunc.clone()
    }
    /// Indicates what “pass” of the pipeline the hit was detected on. This is probably only useful for testing and
    /// debugging. Non-truncated hits are found on the first pass, truncated hits are found on successive passes.
    pub fn pass(&self) -> i32 {
        self.pass
    }
    /// Fraction of G and C nucleotides in the hit.
    pub fn gc(&self) -> f32 {
        self.gc
    }
    /// The biased-composition correction: the bit score difference contributed by the null3 model for CM hits, or
    /// the null2 model for HMM hits. High bias scores may be a red flag for a false positive. It is difficult to correct for all
    /// possible ways in which a nonrandom but nonhomologous biological sequences can appear to be similar, such as
    /// short-period tandem repeats, so there are cases where the bias correction is not strong enough (creating false
    /// positives).
    pub fn bias(&self) -> f32 {
        self.bias
    }
    /// The score (in bits) for this target/query comparison. It includes the biased-composition correction (the
    /// “null3” model for CM hits, or the “null2” model for HMM hits).
    pub fn score(&self) -> f32 {
        self.score
    }
    /// The expectation value (statistical significance) of the target. This is a per query E-value; i.e. calcu-
    /// lated as the expected number of false positives achieving this comparison’s score for a single query against the
    /// search space Z. For cmsearch Z is defined as the total number of nucleotides in the target dataset multiplied
    /// by 2 because both strands are searched. For cmscan Z is the total number of nucleotides in the query sequence
    /// multiplied by 2 because both strands are searched and multiplied by the number of models in the target database.
    /// If you search with multiple queries and if you want to control the overall false positive rate of that search rather
    /// than the false positive rate per query, you will want to multiply this per-query E-value by how many queries you’re
    /// doing.
    pub fn e_value(&self) -> f32 {
        self.e_value
    }
    /// Indicates whether or not this hit achieves the inclusion threshold: ’!’ if it does, ’?’ if it does not (and rather
    /// only achieves the reporting threshold). By default, the inclusion threshold is an E-value of 0.01 and the reporting
    /// threshold is an E-value of 10.0, but these can be changed with command line options as described in the manual
    /// pages.
    pub fn inc(&self) -> char {
        self.inc
    }
    /// The remainder of the line is the target’s description line, as free text.
    pub fn description(&self) -> String {
        self.description.clone()
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
