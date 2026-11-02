//! The _.sufr_ file layout

use crate::{
    io::write_all_at,
    types::{SeedMask, SuffixSortType, OUTFILE_VERSION},
};
use anyhow::{anyhow, bail, Result};
use std::{
    fs::File,
    io::{Read, Seek, SeekFrom},
    mem,
};

// the file is laid out as below, where u64 fields are
// little-endian, and the sequence starts, suffix array
// and LCP array hold `int_size` (4 or 8) bytes per value
// in native byte order
//
//               offset  field
//                    0  version (u8)
//                    1  is_dna (u8, 0/1)
//                    2  allow_ambiguity (u8, 0/1)
//                    3  ignore_softmask (u8, 0/1)
//                    4  text_len (u64)
//                   12  text_pos (u64)
//                   20  suffix_array_pos (u64)
//                   28  lcp_pos (u64, 0 when there is no LCP array)
//                   36  len_suffixes (u64)
//                   44  max_query_len (u64, 0 for a seed mask)
//                   52  num_sequences (u64)
//                   60  sequence starts (num_sequences values)
//                    .  seed mask length (u64)
//                    .  seed mask (u8 per position, 0/1)
//             text_pos  text (text_len bytes)
//     suffix_array_pos  suffix array (len_suffixes values)
//              lcp_pos  LCP array (len_suffixes values), optional
//                    .  sequence names (bincode Vec<String>)
const VERSION: usize = 0;
const IS_DNA: usize = 1;
const ALLOW_AMBIGUITY: usize = 2;
const IGNORE_SOFTMASK: usize = 3;
const TEXT_LEN: usize = 4;
const TEXT_POS: usize = 12;
const SUFFIX_ARRAY_POS: usize = 20;
const LCP_POS: usize = 28;
const LEN_SUFFIXES: usize = 36;
const MAX_QUERY_LEN: usize = 44;
const NUM_SEQUENCES: usize = 52;
const SEQUENCE_STARTS: usize = 60;

// --------------------------------------------------
/// The header of a _.sufr_ file.
#[derive(Debug, PartialEq)]
pub(crate) struct Header {
    pub version: u8,
    pub is_dna: bool,
    pub allow_ambiguity: bool,
    pub ignore_softmask: bool,
    pub text_len: usize,
    pub len_suffixes: usize,
    pub sort_type: SuffixSortType,
    pub sequence_starts: Vec<usize>,
    pub sequence_names: Vec<String>,

    /// The bytes per sequence start, suffix, and LCP value: 4 or 8.
    pub int_size: usize,

    pub text_pos: usize,
    pub suffix_array_pos: usize,

    /// The byte position of the LCP array, if the file has one.
    pub lcp_pos: Option<usize>,
}

impl Header {
    /// Read the header of `filename`, whose arrays hold `int_size`
    /// bytes per value.
    pub fn read(filename: &str, int_size: usize) -> Result<Header> {
        let mut file = File::open(filename).map_err(|e| anyhow!("{filename}: {e}"))?;

        let mut fixed = [0u8; SEQUENCE_STARTS];
        file.read_exact(&mut fixed)?;
        let u64_at =
            |pos: usize| u64::from_le_bytes(fixed[pos..pos + 8].try_into().unwrap());

        let text_len = u64_at(TEXT_LEN);
        if mem::size_of::<usize>() == 4 && text_len > usize::MAX as u64 {
            bail!("text_len exceeds native usize");
        }

        let num_sequences = u64_at(NUM_SEQUENCES) as usize;
        let mut starts = vec![0u8; num_sequences * int_size];
        file.read_exact(&mut starts)?;
        let sequence_starts = starts
            .chunks_exact(int_size)
            .map(|bytes| match int_size {
                4 => u32::from_ne_bytes(bytes.try_into().unwrap()) as usize,
                _ => u64::from_ne_bytes(bytes.try_into().unwrap()) as usize,
            })
            .collect();

        let mut buffer = [0u8; 8];
        file.read_exact(&mut buffer)?;
        let mut seed_mask = vec![0u8; u64::from_le_bytes(buffer) as usize];
        file.read_exact(&mut seed_mask)?;
        let sort_type = if seed_mask.is_empty() {
            SuffixSortType::MaxQueryLen(u64_at(MAX_QUERY_LEN) as usize)
        } else {
            SuffixSortType::Mask(SeedMask::from_bytes(&seed_mask)?)
        };

        let lcp_pos = u64_at(LCP_POS) as usize;
        let mut header = Header {
            version: fixed[VERSION],
            is_dna: fixed[IS_DNA] == 1,
            allow_ambiguity: fixed[ALLOW_AMBIGUITY] == 1,
            ignore_softmask: fixed[IGNORE_SOFTMASK] == 1,
            text_len: text_len as usize,
            len_suffixes: u64_at(LEN_SUFFIXES) as usize,
            sort_type,
            sequence_starts,
            sequence_names: vec![],
            int_size,
            text_pos: u64_at(TEXT_POS) as usize,
            suffix_array_pos: u64_at(SUFFIX_ARRAY_POS) as usize,
            lcp_pos: (lcp_pos > 0).then_some(lcp_pos),
        };

        file.seek(SeekFrom::Start(header.names_pos() as u64))?;
        let mut names = vec![];
        file.read_to_end(&mut names)?;
        header.sequence_names = bincode::deserialize(&names)?;

        Ok(header)
    }

    /// Set the text, suffix array, and LCP positions from the
    /// header's contents, placing an LCP array when `write_lcp`.
    pub fn set_positions(&mut self, write_lcp: bool) {
        let seed_mask_len = match &self.sort_type {
            SuffixSortType::Mask(seed_mask) => seed_mask.bytes.len(),
            SuffixSortType::MaxQueryLen(_) => 0,
        };
        self.text_pos = SEQUENCE_STARTS
            + self.sequence_starts.len() * self.int_size
            + mem::size_of::<u64>()
            + seed_mask_len;
        self.suffix_array_pos = self.text_pos + self.text_len;
        let suffix_array_end =
            self.suffix_array_pos + self.len_suffixes * self.int_size;
        self.lcp_pos = write_lcp.then_some(suffix_array_end);
    }

    /// The byte position of the sequence names, after the arrays.
    pub fn names_pos(&self) -> usize {
        let num_arrays = if self.lcp_pos.is_some() { 2 } else { 1 };
        self.suffix_array_pos + num_arrays * self.len_suffixes * self.int_size
    }

    /// Write the header and the sequence names into `file`, leaving
    /// the text and arrays to the caller. Returns the file's length.
    pub fn write(&self, file: &File) -> Result<usize> {
        let (max_query_len, seed_mask): (usize, &[u8]) = match &self.sort_type {
            SuffixSortType::MaxQueryLen(max_query_len) => (*max_query_len, &[]),
            SuffixSortType::Mask(seed_mask) => (0, &seed_mask.bytes),
        };

        let mut bytes = vec![
            self.version,
            self.is_dna as u8,
            self.allow_ambiguity as u8,
            self.ignore_softmask as u8,
        ];
        bytes.extend((self.text_len as u64).to_le_bytes());
        bytes.extend((self.text_pos as u64).to_le_bytes());
        bytes.extend((self.suffix_array_pos as u64).to_le_bytes());
        bytes.extend((self.lcp_pos.unwrap_or(0) as u64).to_le_bytes());
        bytes.extend((self.len_suffixes as u64).to_le_bytes());
        bytes.extend((max_query_len as u64).to_le_bytes());
        bytes.extend((self.sequence_starts.len() as u64).to_le_bytes());
        for &start in &self.sequence_starts {
            match self.int_size {
                4 => bytes.extend((start as u32).to_ne_bytes()),
                _ => bytes.extend((start as u64).to_ne_bytes()),
            }
        }
        bytes.extend((seed_mask.len() as u64).to_le_bytes());
        bytes.extend_from_slice(seed_mask);
        debug_assert_eq!(bytes.len(), self.text_pos);
        write_all_at(file, &bytes, 0)?;

        let names = bincode::serialize(&self.sequence_names)?;
        write_all_at(file, &names, self.names_pos() as u64)?;

        Ok(self.names_pos() + names.len())
    }
}

// --------------------------------------------------
/// Read the length of the text in a _.sufr_ file, which decides
/// whether its arrays hold `u32` or `u64` values.
pub fn read_text_length(filename: &str) -> Result<usize> {
    let mut file = File::open(filename).map_err(|e| anyhow!("{filename}: {e}"))?;

    let mut flags = [0u8; TEXT_LEN];
    file.read_exact(&mut flags)?;
    let version = flags[VERSION];
    if version != OUTFILE_VERSION {
        bail!("Unknown sufr version {version}");
    }

    let mut buffer = [0u8; 8];
    file.read_exact(&mut buffer)?;
    let text_len = u64::from_le_bytes(buffer);
    if mem::size_of::<usize>() == 4 && text_len > usize::MAX as u64 {
        bail!("text_len exceeds native usize");
    }
    Ok(text_len as usize)
}

// --------------------------------------------------
#[cfg(test)]
mod tests {
    use super::{read_text_length, Header};
    use crate::types::{SeedMask, SuffixSortType, OUTFILE_VERSION};
    use anyhow::Result;
    use pretty_assertions::assert_eq;
    use tempfile::NamedTempFile;

    #[test]
    fn test_read_text_length() -> Result<()> {
        let sufr_file = "../data/inputs/2.sufr";
        let res = read_text_length(sufr_file);
        assert!(res.is_ok());
        let len = res.unwrap();
        assert_eq!(len, 18);
        Ok(())
    }

    #[test]
    fn test_header_round_trip() -> Result<()> {
        let cases = [
            (4, SuffixSortType::Mask(SeedMask::new("1101")?), true),
            (8, SuffixSortType::MaxQueryLen(5), false),
        ];
        for (int_size, sort_type, write_lcp) in cases {
            let mut header = Header {
                version: OUTFILE_VERSION,
                is_dna: true,
                allow_ambiguity: false,
                ignore_softmask: true,
                text_len: 20,
                len_suffixes: 17,
                sort_type,
                sequence_starts: vec![0, 7, 13],
                sequence_names: vec!["a".to_string(), "b".into(), "c".into()],
                int_size,
                text_pos: 0,
                suffix_array_pos: 0,
                lcp_pos: None,
            };
            header.set_positions(write_lcp);

            let outfile = NamedTempFile::new()?;
            let len = header.write(outfile.as_file())?;
            assert_eq!(len as u64, outfile.as_file().metadata()?.len());

            let path = outfile.path().to_string_lossy().to_string();
            assert_eq!(Header::read(&path, int_size)?, header);
        }
        Ok(())
    }
}
