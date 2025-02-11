use std::{
    fs::File,
    io::{self, BufWriter, Write},
};

use crate::{record::Header, Meta, Record};

/// Writer for HMMER-compatible files
pub struct Writer<W: Write> {
    writer: BufWriter<W>,
}

impl<W: Write> Writer<W> {
    /// Creates a new Writer instance
    pub fn new(writer: W) -> Self {
        Self {
            writer: BufWriter::new(writer),
        }
    }

    /// Writes a meta line in HMMER-compatible format
    pub fn write_meta(&mut self, meta: Meta) -> io::Result<()> {
        write!(self.writer, "{}", meta)
    }

    /// Writes a program line in HMMER-compatible format
    pub fn write_header(&mut self, header: Header) -> io::Result<()> {
        writeln!(self.writer, "{}", header)
    }

    /// Writes a single record in HMMER-compatible format
    pub fn write_record(&mut self, record: &Record) -> io::Result<()> {
        match record {
            Record::Dna(record) => writeln!(self.writer, "{}", record),
            Record::Protein(record) => writeln!(self.writer, "{}", record),
            Record::CM(cmrecord) => writeln!(self.writer, "{}", cmrecord),
        }
    }

    /// Flushes the buffer to ensure all data is written
    pub fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }

    /// Consumes the Writer and returns the inner writer
    pub fn into_inner(self) -> Result<W, io::IntoInnerError<BufWriter<W>>> {
        self.writer.into_inner()
    }
}

impl Writer<File> {
    /// Convenience function to create a Writer that writes to a file
    pub fn to_file(path: &str) -> io::Result<Self> {
        let file = File::create(path)?;
        Ok(Self::new(file))
    }
}
