use crate::{
    record::{CMRecord, Header, Meta, Program, Record, Strand},
    DNARecord, Error, ErrorKind, ProteinRecord, Result,
};

use std::{
    collections::VecDeque,
    fs::File,
    io::{self, BufRead},
    path::{Path, PathBuf},
    str::FromStr,
};

/// A reader over the header of a HMM tblout file.
pub struct HeaderReader<R> {
    rdr: io::BufReader<R>,
    line: u64,
}

impl<R: io::Read> HeaderReader<R> {
    pub fn new(rdr: R) -> HeaderReader<R> {
        HeaderReader {
            rdr: io::BufReader::new(rdr),
            line: 0,
        }
    }

    fn read_header(&mut self) -> Result<Header> {
        let mut header = Header::default();

        let mut line = String::new();
        loop {
            line.clear();
            match self.rdr.read_line(&mut line) {
                Ok(0) => break,
                Ok(_) => {
                    self.line += 1;
                    let mut hash_counter = 1;

                    if line.starts_with('#') {
                        // probably more checks here than necessary.
                        if line.contains("full sequence") {
                            header.set_protein_only(line.clone());
                            continue;
                        }

                        if line.contains("target name") {
                            header.set_columns(line.clone());
                            continue;
                        }

                        // Check here that the header only contains dashes
                        fn remove_whitespace(s: &str) -> String {
                            s.chars().filter(|c| !c.is_whitespace()).collect()
                        }

                        fn is_single_char(s: &str) -> bool {
                            if let Some(first) = s.chars().next() {
                                s.chars().all(|c| c == first)
                            } else {
                                false // Empty string case
                            }
                        }

                        let dashes_line_prefix_removed = line.strip_prefix('#').unwrap();
                        let trimmed_line = remove_whitespace(dashes_line_prefix_removed);

                        if is_single_char(&trimmed_line) {
                            header.set_dashes(line.clone());
                            continue;
                        }

                        hash_counter += 1;
                    }

                    if hash_counter > 3 {
                        break;
                    }
                }
                Err(e) => return Err(Error::new(ErrorKind::Io(e))),
            }
        }

        Ok(header)
    }
}

/// A reader over the metadata of a HMM tblout file.
pub struct MetaReader<R> {
    rdr: io::BufReader<R>,
    line: u64,
}

impl<R: io::Read> MetaReader<R> {
    /// Constructor for `MetaReader`.
    pub fn new(rdr: R) -> MetaReader<R> {
        MetaReader {
            rdr: io::BufReader::new(rdr),
            line: 0,
        }
    }
    /// Read the metadata from the input reader.
    fn read_meta(&mut self) -> Result<Meta> {
        // read the metadata into the meta struct
        // we skip the first three #'s that we come across
        // and the fourth should be where the metadata starts
        let mut metadata = Meta::default();

        let mut line = String::new();
        let mut hash_counter = 0;
        loop {
            line.clear();
            match self.rdr.read_line(&mut line) {
                Ok(0) => break,
                Ok(_) => {
                    self.line += 1;

                    // increment the hash counter
                    if line.starts_with('#') {
                        hash_counter += 1;
                    }

                    if !line.starts_with('#') {
                        continue;
                    }

                    // once we hit the fourth hash we can start reading
                    if hash_counter >= 4 {
                        // match on the starting substring of the line
                        let mut split_line = line
                            // FIXME: if the line contains other instances of ": ",
                            // some of this will fail to parse correctly
                            .split(": ")
                            .collect::<Vec<&str>>()
                            .into_iter()
                            .map(|e| e.trim())
                            .collect::<VecDeque<&str>>();

                        let first = split_line.pop_front().unwrap();
                        let rest = split_line.into_iter().collect::<Vec<&str>>().join(" ");

                        match first {
                            "# Program" => metadata.set_program(Program::from_str(&rest).unwrap()),
                            "# Version" => metadata.set_version(rest.to_string()),
                            "# Pipeline mode" => metadata.set_pipeline_mode(rest.to_string()),
                            "# Query file" => {
                                metadata.set_query_file(PathBuf::from(rest.to_string()))
                            }
                            "# Target file" => {
                                metadata.set_target_file(PathBuf::from(rest.to_string()))
                            }
                            "# Option settings" => metadata.set_options(rest.to_string()),
                            "# Current dir" => {
                                metadata.set_current_dir(PathBuf::from(rest.to_string()))
                            }
                            "# Date" => metadata.set_date(rest.to_string()),
                            _ => (), // probably make this an error
                        }
                    }
                }

                Err(e) => return Err(Error::new(ErrorKind::Io(e))),
            }
        }

        Ok(metadata)
    }
}

/// A reader over the records of a HMM tblout file.
pub struct Reader<R> {
    /// The underlying reader.
    rdr: io::BufReader<R>,
    /// The current line number.
    line: u64,
    /// The metadata from the first pass.
    meta: Meta,
    /// The header from the first pass.
    header: Header,
}

impl Reader<File> {
    /// Construct a new reader from a file path.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Reader<File>> {
        let mut headerreader = HeaderReader::new(File::open(path.as_ref())?);
        let header = headerreader.read_header()?;

        let mut metareader = MetaReader::new(File::open(path.as_ref())?);
        let meta = metareader.read_meta()?;

        Ok(Reader::new(File::open(path)?, header, meta))
    }

    /// Construct a new reader from anything that implements `io::Read`
    /// and clone.
    pub fn from_reader<R: io::Read + Clone>(rdr: R) -> Result<Reader<R>> {
        let mut headerreader = HeaderReader::new(rdr.clone());
        let header = headerreader.read_header()?;

        let mut metareader = MetaReader::new(rdr.clone());
        let meta = metareader.read_meta()?;

        Ok(Reader::new(rdr, header, meta))
    }
}

/// An iterator over the records of a refer file.
pub enum RecordsIter<'a, R> {
    /// An iterator over DNA records.
    Dna(DNARecordsIter<'a, R>),
    /// An iterator over protein records.
    Protein(ProteinRecordsIter<'a, R>),
    /// An iterator over CM records.
    CM(CMRecordsIter<'a, R>),
}

impl<'r, R: io::Read> Iterator for RecordsIter<'r, R> {
    type Item = Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            RecordsIter::Dna(e) => e.next().map(|rec| rec.map(Record::Dna)),
            RecordsIter::Protein(e) => e.next().map(|rec| rec.map(Record::Protein)),
            RecordsIter::CM(e) => e.next().map(|rec| rec.map(Record::CM)),
        }
    }
}

/// A borrowed iterator over the records of a refer file.
pub struct DNARecordsIter<'r, R: 'r> {
    /// The underlying reader
    rdr: &'r mut Reader<R>,
}

impl<'r, R: io::Read> Iterator for DNARecordsIter<'r, R> {
    type Item = Result<DNARecord>;

    fn next(&mut self) -> Option<Result<DNARecord>> {
        match self.rdr.read_dna_record() {
            Ok(Some(r)) => {
                self.rdr.line += 1;
                Some(Ok(r))
            }
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

/// A borrowed iterator over the records of a refer file.
pub struct ProteinRecordsIter<'r, R: 'r> {
    /// The underlying reader
    rdr: &'r mut Reader<R>,
}

impl<'r, R: io::Read> Iterator for ProteinRecordsIter<'r, R> {
    type Item = Result<ProteinRecord>;

    fn next(&mut self) -> Option<Result<ProteinRecord>> {
        match self.rdr.read_protein_record() {
            Ok(Some(r)) => {
                self.rdr.line += 1;
                Some(Ok(r))
            }
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

/// A borrowed iterator over the records of a refer file.
pub struct CMRecordsIter<'r, R: 'r> {
    /// The underlying reader
    rdr: &'r mut Reader<R>,
}

impl<'r, R: io::Read> Iterator for CMRecordsIter<'r, R> {
    type Item = Result<CMRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.rdr.read_cm_record() {
            Ok(Some(r)) => {
                self.rdr.line += 1;
                Some(Ok(r))
            }
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

impl<'r, R: io::Read> RecordsIter<'r, R> {
    fn new(rdr: &'r mut Reader<R>, program: Program) -> RecordsIter<'r, R> {
        match program {
            Program::Nhmmer | Program::Nhmmscan => RecordsIter::Dna(DNARecordsIter { rdr }),
            Program::Jackhmmer | Program::Hmmscan | Program::Hmmsearch | Program::Phmmer => {
                RecordsIter::Protein(ProteinRecordsIter { rdr })
            }
            Program::Cmsearch | Program::Cmscan => RecordsIter::CM(CMRecordsIter { rdr }),
            Program::None => unreachable!(),
        }
    }
    /// Return a reference to the underlying reader.
    pub fn reader(&self) -> &Reader<R> {
        match self {
            RecordsIter::Dna(r) => r.rdr,
            RecordsIter::Protein(r) => r.rdr,
            RecordsIter::CM(r) => r.rdr,
        }
    }

    /// Return a mutable reference to the underlying reader.
    pub fn reader_mut(&mut self) -> &mut Reader<R> {
        match self {
            RecordsIter::Dna(r) => r.rdr,
            RecordsIter::Protein(r) => r.rdr,
            RecordsIter::CM(r) => r.rdr,
        }
    }
}

/// An owned iterator over the records of a refer file.
pub enum RecordsIntoIter<R> {
    /// An iterator over DNA records.
    Dna(DNARecordsIntoIter<R>),
    /// An iterator over protein records.
    Protein(ProteinRecordsIntoIter<R>),
    /// An iterator over CM records.
    CM(CMRecordsIntoIter<R>),
}

impl<R: io::Read> Iterator for RecordsIntoIter<R> {
    type Item = Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            RecordsIntoIter::Dna(e) => e.next().map(|rec| rec.map(Record::Dna)),
            RecordsIntoIter::Protein(e) => e.next().map(|rec| rec.map(Record::Protein)),
            RecordsIntoIter::CM(e) => e.next().map(|rec| rec.map(Record::CM)),
        }
    }
}
/// An owned iterator over the records of a DNA tblout file.
pub struct DNARecordsIntoIter<R> {
    /// The underlying reader.
    rdr: Reader<R>,
}

impl<R: io::Read> RecordsIntoIter<R> {
    fn new(rdr: Reader<R>, program: Program) -> RecordsIntoIter<R> {
        match program {
            Program::Nhmmer | Program::Nhmmscan => RecordsIntoIter::Dna(DNARecordsIntoIter { rdr }),
            Program::Jackhmmer | Program::Hmmscan | Program::Hmmsearch | Program::Phmmer => {
                RecordsIntoIter::Protein(ProteinRecordsIntoIter { rdr })
            }
            Program::Cmsearch | Program::Cmscan => RecordsIntoIter::CM(CMRecordsIntoIter { rdr }),
            Program::None => unreachable!(),
        }
    }
    /// Return a reference to the underlying reader.
    pub fn reader(&self) -> &Reader<R> {
        match self {
            RecordsIntoIter::Dna(r) => &r.rdr,
            RecordsIntoIter::Protein(r) => &r.rdr,
            RecordsIntoIter::CM(r) => &r.rdr,
        }
    }

    /// Return a mutable reference to the underlying reader.
    pub fn reader_mut(&mut self) -> &mut Reader<R> {
        match self {
            RecordsIntoIter::Dna(r) => &mut r.rdr,
            RecordsIntoIter::Protein(r) => &mut r.rdr,
            RecordsIntoIter::CM(r) => &mut r.rdr,
        }
    }

    /// Drop this iterator and return the underlying reader.
    pub fn into_reader(self) -> Reader<R> {
        match self {
            RecordsIntoIter::Dna(r) => r.rdr,
            RecordsIntoIter::Protein(r) => r.rdr,
            RecordsIntoIter::CM(r) => r.rdr,
        }
    }
}

impl<R: io::Read> Iterator for DNARecordsIntoIter<R> {
    type Item = Result<DNARecord>;

    fn next(&mut self) -> Option<Result<DNARecord>> {
        match self.rdr.read_dna_record() {
            Ok(Some(r)) => {
                self.rdr.line += 1;
                Some(Ok(r))
            }
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

/// An owned iterator over the records of a protein tblout file.
pub struct ProteinRecordsIntoIter<R> {
    /// The underlying reader.
    rdr: Reader<R>,
}

impl<R: io::Read> Iterator for ProteinRecordsIntoIter<R> {
    type Item = Result<ProteinRecord>;

    fn next(&mut self) -> Option<Result<ProteinRecord>> {
        match self.rdr.read_protein_record() {
            Ok(Some(r)) => {
                self.rdr.line += 1;
                Some(Ok(r))
            }
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

/// An owned iterator over the records of a CM tblout file.
pub struct CMRecordsIntoIter<R> {
    /// The underlying reader.
    rdr: Reader<R>,
}

impl<R: io::Read> Iterator for CMRecordsIntoIter<R> {
    type Item = Result<CMRecord>;

    fn next(&mut self) -> Option<Result<CMRecord>> {
        match self.rdr.read_cm_record() {
            Ok(Some(r)) => {
                self.rdr.line += 1;
                Some(Ok(r))
            }
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

impl<R: io::Read> Reader<R> {
    /// Construct a new reader from a reader and metadata.
    pub fn new(rdr: R, header: Header, meta: Meta) -> Reader<R> {
        Reader {
            rdr: io::BufReader::new(rdr),
            line: 0,
            header,
            meta,
        }
    }

    /// Return the header from the first pass.
    pub fn header(&self) -> &Header {
        &self.header
    }

    /// Return the metadata from the first pass.
    pub fn meta(&self) -> &Meta {
        &self.meta
    }

    /// A borrowed iterator over the records of a refer file.
    pub fn records(&mut self) -> RecordsIter<R> {
        RecordsIter::new(self, self.meta.program())
    }

    /// An owned iterator over the records of a refer file.
    pub fn into_records(self) -> RecordsIntoIter<R> {
        let program = self.meta.program();
        RecordsIntoIter::new(self, program)
    }

    /// Read a single record from an input reader.
    fn read_dna_record(&mut self) -> Result<Option<DNARecord>> {
        // for this function, we read a single line and parse
        // on whitespace, returning a record. We skip lines
        // starting with a comment character '#'.
        let mut line = String::new();
        let col_sizes = self.header.calculate_dashes();

        loop {
            line.clear();
            match self.rdr.read_line(&mut line) {
                Ok(0) => return Ok(None),
                Ok(_) => {
                    self.line += 1;
                    if line.starts_with('#') {
                        continue;
                    }
                    let l_vec = line.split_whitespace().collect::<Vec<&str>>();

                    let target_name = l_vec[0].to_string();
                    let target_accession = l_vec[1].to_string();
                    let query_name = l_vec[2].to_string();
                    let query_accession = l_vec[3].to_string();
                    let hmm_from = l_vec[4].parse::<i32>()?;
                    let hmm_to = l_vec[5].parse::<i32>()?;
                    let ali_from = l_vec[6].parse::<i32>()?;
                    let ali_to = l_vec[7].parse::<i32>()?;
                    let env_from = l_vec[8].parse::<i32>()?;
                    let env_to = l_vec[9].parse::<i32>()?;
                    let sq_len = l_vec[10].parse::<i32>()?;
                    let strand = l_vec[11].parse::<Strand>()?;
                    let e_value = l_vec[12].parse::<f32>()?;
                    let score = l_vec[13].parse::<f32>()?;
                    let bias = l_vec[14].parse::<f32>()?;

                    let description = l_vec[15..].join(" ");

                    let record = DNARecord::new(
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
                    );

                    return Ok(Some(record));
                }
                Err(e) => return Err(Error::new(ErrorKind::Io(e))),
            }
        }
    }

    /// Read a single record from an input reader as a protein record.
    fn read_protein_record(&mut self) -> Result<Option<ProteinRecord>> {
        // for this function, we read a single line and parse
        // on whitespace, returning a record. We skip lines
        // starting with a comment character '#'.
        let mut line = String::new();
        let col_sizes = self.header.calculate_dashes();
        loop {
            line.clear();
            match self.rdr.read_line(&mut line) {
                Ok(0) => return Ok(None),
                Ok(_) => {
                    self.line += 1;
                    if line.starts_with('#') {
                        continue;
                    }
                    let l_vec = line.split_whitespace().collect::<Vec<&str>>();
                    let target_name = l_vec[0].to_string();
                    let target_accession = l_vec[1].to_string();
                    let query_name = l_vec[2].to_string();
                    let query_accession = l_vec[3].to_string();
                    let e_value_full = l_vec[4].parse::<f32>()?;
                    let score_full = l_vec[5].parse::<f32>()?;
                    let bias_full = l_vec[6].parse::<f32>()?;
                    let e_value_best = l_vec[7].parse::<f32>()?;
                    let score_best = l_vec[8].parse::<f32>()?;
                    let bias_best = l_vec[9].parse::<f32>()?;
                    let exp = l_vec[10].parse::<f32>()?;
                    let reg = l_vec[11].parse::<i32>()?;
                    let clu = l_vec[12].parse::<i32>()?;
                    let ov = l_vec[13].parse::<i32>()?;
                    let env = l_vec[14].parse::<i32>()?;
                    let dom = l_vec[15].parse::<i32>()?;
                    let rep = l_vec[16].parse::<i32>()?;
                    let inc = l_vec[17].parse::<i32>()?;

                    // description is all the remainder of the line as free text
                    let description = l_vec[18..].join(" ");

                    let record = ProteinRecord::new(
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
                    );

                    return Ok(Some(record));
                }
                Err(e) => return Err(Error::new(ErrorKind::Io(e))),
            }
        }
    }

    fn read_cm_record(&mut self) -> Result<Option<CMRecord>> {
        let mut line = String::new();
        let col_sizes = self.header.calculate_dashes();
        loop {
            line.clear();
            match self.rdr.read_line(&mut line) {
                Ok(0) => return Ok(None),
                Ok(_) => {
                    self.line += 1;
                    if line.starts_with('#') {
                        continue;
                    }
                    let l_vec = line.split_whitespace().collect::<Vec<&str>>();
                    let target_name = l_vec[0].to_string();
                    let target_accession = l_vec[1].to_string();
                    let query_name = l_vec[2].to_string();
                    let query_accession = l_vec[3].to_string();

                    let mdl = l_vec[4].to_string();
                    let mdl_from = l_vec[5].parse::<i32>()?;
                    let mdl_to = l_vec[6].parse::<i32>()?;

                    let seq_from = l_vec[7].parse::<i32>()?;
                    let seq_to = l_vec[8].parse::<i32>()?;

                    let strand = l_vec[9].parse::<Strand>()?;

                    let trunc = l_vec[10].to_string();
                    let pass = l_vec[11].parse::<i32>()?;

                    let gc = l_vec[12].parse::<f32>()?;

                    let bias = l_vec[13].parse::<f32>()?;
                    let score = l_vec[14].parse::<f32>()?;

                    let e_value = l_vec[15].parse::<f32>()?;

                    let inc = l_vec[16].parse::<char>()?;

                    // description is all the remainder of the line as free text
                    let description = l_vec[17..].join(" ");

                    let record = CMRecord::new(
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
                    );

                    return Ok(Some(record));
                }
                Err(e) => return Err(Error::new(ErrorKind::Io(e))),
            }
        }
    }
}
