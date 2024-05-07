use crate::{
    record::{Meta, Program, Strand},
    Error, ErrorKind, Record, Result,
};

use std::{
    collections::VecDeque,
    fs::File,
    io::{self, BufRead},
    path::{Path, PathBuf},
    str::FromStr,
};

pub struct MetaReader<R> {
    rdr: io::BufReader<R>,
    line: u64,
}

impl<R: io::Read> MetaReader<R> {
    pub fn new(rdr: R) -> MetaReader<R> {
        MetaReader {
            rdr: io::BufReader::new(rdr),
            line: 0,
        }
    }
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
                            .split(':')
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

pub struct Reader<R> {
    rdr: io::BufReader<R>,
    line: u64,
    meta: Meta,
}

impl Reader<File> {
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Reader<File>> {
        let mut metareader = MetaReader::new(File::open(path.as_ref())?);
        let meta = metareader.read_meta()?;

        Ok(Reader::new(File::open(path)?, meta))
    }
    pub fn from_reader<R: io::Read + Clone>(rdr: R) -> Result<Reader<R>> {
        let mut metareader = MetaReader::new(rdr.clone());
        let meta = metareader.read_meta()?;

        Ok(Reader::new(rdr, meta))
    }
}

impl<R: io::Read> Reader<R> {
    pub fn new(rdr: R, meta: Meta) -> Reader<R> {
        Reader {
            rdr: io::BufReader::new(rdr),
            line: 0,
            meta,
        }
    }

    /// Return the metadata from the first pass
    pub fn meta(&self) -> &Meta {
        &self.meta
    }

    /// A borrowed iterator over the records of a refer file.
    pub fn records(&mut self) -> RecordsIter<R> {
        RecordsIter::new(self)
    }

    /// An owned iterator over the records of a refer file.
    pub fn into_records(self) -> RecordsIntoIter<R> {
        RecordsIntoIter::new(self)
    }

    /// Read a single record from an input reader.
    fn read_record(&mut self) -> Result<Option<Record>> {
        // reset the line counter
        self.line = 0;
        // for this function, we read a single line and parse
        // on whitespace, returning a record. We skip lines
        // starting with a comment character '#'.
        let mut line = String::new();
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
                    let query_name = l_vec[2].to_string();
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
                    // note we omit description column

                    let record = Record::new(
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
                    );

                    return Ok(Some(record));
                }
                Err(e) => return Err(Error::new(ErrorKind::Io(e))),
            }
        }
    }
}

/// A borrowed iterator over the records of a refer file.
pub struct RecordsIter<'r, R: 'r> {
    /// The underlying reader
    rdr: &'r mut Reader<R>,
}

impl<'r, R: io::Read> RecordsIter<'r, R> {
    fn new(rdr: &'r mut Reader<R>) -> RecordsIter<'r, R> {
        RecordsIter { rdr }
    }
    /// Return a reference to the underlying reader.
    pub fn reader(&self) -> &Reader<R> {
        self.rdr
    }

    /// Return a mutable reference to the underlying reader.
    pub fn reader_mut(&mut self) -> &mut Reader<R> {
        self.rdr
    }
}

impl<'r, R: io::Read> Iterator for RecordsIter<'r, R> {
    type Item = Result<Record>;

    fn next(&mut self) -> Option<Result<Record>> {
        match self.rdr.meta.program() {
            Program::None => Some(Err(Error::new(ErrorKind::Meta("program not read".into())))),
            Program::Nhmmer => match self.rdr.read_record() {
                Ok(Some(r)) => {
                    self.rdr.line += 1;
                    Some(Ok(r))
                }
                Ok(None) => None,
                Err(e) => Some(Err(e)),
            },
        }
    }
}

/// An owned iterator over the records of a refer file.
pub struct RecordsIntoIter<R> {
    /// The underlying reader.
    rdr: Reader<R>,
}

impl<R: io::Read> RecordsIntoIter<R> {
    fn new(rdr: Reader<R>) -> RecordsIntoIter<R> {
        RecordsIntoIter { rdr }
    }
    /// Return a reference to the underlying reader.
    pub fn reader(&self) -> &Reader<R> {
        &self.rdr
    }

    /// Return a mutable reference to the underlying reader.
    pub fn reader_mut(&mut self) -> &mut Reader<R> {
        &mut self.rdr
    }

    /// Drop this iterator and return the underlying reader.
    pub fn into_reader(self) -> Reader<R> {
        self.rdr
    }
}

impl<R: io::Read> Iterator for RecordsIntoIter<R> {
    type Item = Result<Record>;

    fn next(&mut self) -> Option<Result<Record>> {
        match self.rdr.meta.program() {
            Program::None => Some(Err(Error::new(ErrorKind::Meta("program not read".into())))),
            Program::Nhmmer => match self.rdr.read_record() {
                Ok(Some(r)) => {
                    self.rdr.line += 1;
                    Some(Ok(r))
                }
                Ok(None) => None,
                Err(e) => Some(Err(e)),
            },
        }
    }
}
