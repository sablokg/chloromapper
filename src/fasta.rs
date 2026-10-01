//! Minimal, dependency-free FASTA I/O.

use std::fmt;
use std::fs;
use std::io::{self, Write};
use std::path::Path;

/*
Gaurav Sablok
gsablok@proton.me
*/

#[derive(Debug, Clone)]
pub struct FastaRecord {
    /// Text after '>' up to the first whitespace.
    pub id: String,
    /// Remainder of the header line after the id (may be empty).
    pub description: String,
    /// Raw sequence bytes, newlines stripped, case preserved.
    pub seq: Vec<u8>,
}

impl FastaRecord {
    pub fn header(&self) -> String {
        if self.description.is_empty() {
            self.id.clone()
        } else {
            format!("{} {}", self.id, self.description)
        }
    }
}

#[derive(Debug)]
pub enum FastaError {
    Io(io::Error),
    Empty,
    NoHeader(usize),
}

impl fmt::Display for FastaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FastaError::Io(e) => write!(f, "I/O error: {e}"),
            FastaError::Empty => write!(f, "input contains no FASTA records"),
            FastaError::NoHeader(line) => {
                write!(f, "sequence data before any '>' header at line {line}")
            }
        }
    }
}

impl std::error::Error for FastaError {}

impl From<io::Error> for FastaError {
    fn from(e: io::Error) -> Self {
        FastaError::Io(e)
    }
}

/// Parse FASTA text into records. Tolerates '\r\n', blank lines, and lowercase bases.
pub fn parse(text: &str) -> Result<Vec<FastaRecord>, FastaError> {
    let mut records = Vec::new();
    let mut cur_id: Option<String> = None;
    let mut cur_desc = String::new();
    let mut cur_seq: Vec<u8> = Vec::new();

    for (lineno, raw_line) in text.lines().enumerate() {
        let line = raw_line.trim_end_matches('\r');
        if line.is_empty() {
            continue;
        }
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(id) = cur_id.take() {
                records.push(FastaRecord {
                    id,
                    description: std::mem::take(&mut cur_desc),
                    seq: std::mem::take(&mut cur_seq),
                });
            }
            let mut parts = rest.splitn(2, char::is_whitespace);
            cur_id = Some(parts.next().unwrap_or("").to_string());
            cur_desc = parts.next().unwrap_or("").trim().to_string();
        } else {
            if cur_id.is_none() {
                return Err(FastaError::NoHeader(lineno + 1));
            }
            cur_seq.extend(line.bytes().filter(|b| !b.is_ascii_whitespace()));
        }
    }
    if let Some(id) = cur_id.take() {
        records.push(FastaRecord {
            id,
            description: cur_desc,
            seq: cur_seq,
        });
    }
    if records.is_empty() {
        return Err(FastaError::Empty);
    }
    Ok(records)
}

pub fn read_file<P: AsRef<Path>>(path: P) -> Result<Vec<FastaRecord>, FastaError> {
    let text = fs::read_to_string(path)?;
    parse(&text)
}

/// Write records as wrapped FASTA (default 70 columns).
pub fn write<W: Write>(mut w: W, records: &[FastaRecord], wrap: usize) -> io::Result<()> {
    let wrap = if wrap == 0 { usize::MAX } else { wrap };
    for r in records {
        writeln!(w, ">{}", r.header())?;
        for chunk in r.seq.chunks(wrap) {
            w.write_all(chunk)?;
            w.write_all(b"\n")?;
        }
    }
    Ok(())
}

pub fn write_string(records: &[FastaRecord], wrap: usize) -> String {
    let mut buf = Vec::new();
    write(&mut buf, records, wrap).expect("writing to Vec<u8> cannot fail");
    String::from_utf8(buf).expect("FASTA output is always valid UTF-8/ASCII")
}
