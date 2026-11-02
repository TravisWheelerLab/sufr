//! The _.sufr_ file layout

use crate::io::write_all_at;
use anyhow::{anyhow, bail, Result};
use regex::Regex;
use std::{
    fmt::{self, Display},
    fs::File,
    io::{Read, Seek, SeekFrom},
    mem,
};

/// The serialization version.
pub const OUTFILE_VERSION: u8 = 6;

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
//                   36  num_suffixes (u64)
//                   44  max_query_len (u64, 0 for a seed mask)
//                   52  num_sequences (u64)
//                   60  sequence starts (num_sequences values)
//                    .  seed mask length (u64)
//                    .  seed mask (u8 per position, 0/1)
//             text_pos  text (text_len bytes)
//     suffix_array_pos  suffix array (num_suffixes values)
//              lcp_pos  LCP array (num_suffixes values), optional
//                    .  sequence names (bincode Vec<String>)
const VERSION: usize = 0;
const IS_DNA: usize = 1;
const ALLOW_AMBIGUITY: usize = 2;
const IGNORE_SOFTMASK: usize = 3;
const TEXT_LEN: usize = 4;
const TEXT_POS: usize = 12;
const SUFFIX_ARRAY_POS: usize = 20;
const LCP_POS: usize = 28;
const NUM_SUFFIXES: usize = 36;
const MAX_QUERY_LEN: usize = 44;
const NUM_SEQUENCES: usize = 52;
const SEQUENCE_STARTS: usize = 60;

/// The header of a _.sufr_ file.
#[derive(Debug, PartialEq)]
pub struct Header {
    /// The serialization version.
    pub version: u8,

    /// Whether the text holds nucleotides.
    pub is_dna: bool,

    /// Whether nucleotides other than A, C, G, and T were allowed.
    pub allow_ambiguity: bool,

    /// Whether softmasked (lowercase) nucleotides were ignored.
    pub ignore_softmask: bool,

    /// The length of the text.
    pub text_len: usize,

    /// The number of suffixes in the suffix array.
    pub num_suffixes: usize,

    /// How the suffixes were sorted.
    pub sort_type: SuffixSortType,

    /// The start of each sequence in the text.
    pub sequence_starts: Vec<usize>,

    /// The name of each sequence.
    pub sequence_names: Vec<String>,

    /// The bytes per sequence start, suffix, and LCP value: 4 or 8.
    pub int_size: usize,

    /// The byte position of the text.
    pub text_pos: usize,

    /// The byte position of the suffix array.
    pub suffix_array_pos: usize,

    /// The byte position of the LCP array, if the file has one.
    pub lcp_pos: Option<usize>,
}

impl Header {
    /// Read the header of `filename`, whose arrays hold `int_size`
    /// bytes per value.
    pub(crate) fn read(filename: &str, int_size: usize) -> Result<Header> {
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
            num_suffixes: u64_at(NUM_SUFFIXES) as usize,
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
    pub(crate) fn set_positions(&mut self, write_lcp: bool) {
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
            self.suffix_array_pos + self.num_suffixes * self.int_size;
        self.lcp_pos = write_lcp.then_some(suffix_array_end);
    }

    /// The byte position of the sequence names, after the arrays.
    pub(crate) fn names_pos(&self) -> usize {
        let num_arrays = if self.lcp_pos.is_some() { 2 } else { 1 };
        self.suffix_array_pos + num_arrays * self.num_suffixes * self.int_size
    }

    /// Write the header and the sequence names into `file`, leaving
    /// the text and arrays to the caller. Returns the file's length.
    pub(crate) fn write(&self, file: &File) -> Result<usize> {
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
        bytes.extend((self.num_suffixes as u64).to_le_bytes());
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

/// Read the length of the text in a _.sufr_ file, which decides
/// whether its arrays hold `u32` or `u64` values.
pub(crate) fn read_text_length(filename: &str) -> Result<usize> {
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

/// How the suffixes in a _.sufr_ file were sorted.
#[derive(Debug, Clone, PartialEq)]
pub enum SuffixSortType {
    /// By their first `n` characters, or by all of them when `n` is 0.
    MaxQueryLen(usize),

    /// By the care positions of a seed mask.
    Mask(SeedMask),
}

/// A seed mask of care (1) and don't-care (0) positions.
#[derive(Debug, PartialEq, Clone)]
pub struct SeedMask {
    /// The mask as a string of 1s and 0s.
    pub mask: String,

    /// The mask as one 1 or 0 byte per position.
    pub bytes: Vec<u8>,

    /// The offsets of the care positions.
    pub positions: Vec<usize>,

    /// The difference between each care position and its index in
    /// `positions`.
    pub differences: Vec<usize>,

    /// The number of care positions.
    pub weight: usize,
}

impl SeedMask {
    /// Parse a mask of 1s and 0s that starts and ends with 1 and holds
    /// at least one 0.
    pub fn new(mask: &str) -> Result<Self> {
        if !Self::is_valid(mask) {
            bail!("Invalid seed mask '{mask}'")
        }

        let bytes = Self::parse(mask);
        let positions = Self::get_positions(&bytes);
        let differences = Self::get_differences(&positions);
        let weight = positions.len();

        Ok(Self {
            mask: mask.to_string(),
            bytes,
            positions,
            differences,
            weight,
        })
    }

    /// The mask stored in a _.sufr_ header as one 1 or 0 byte per
    /// position.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self> {
        let mask: String = bytes
            .iter()
            .filter_map(|b| match b {
                0 => Some('0'),
                1 => Some('1'),
                _ => None,
            })
            .collect();
        if mask.is_empty() {
            bail!("Bytes must be 1 or 0");
        }
        let positions = Self::get_positions(bytes);
        let differences = Self::get_differences(&positions);
        let weight = positions.len();

        Ok(Self {
            mask,
            bytes: bytes.to_vec(),
            positions,
            differences,
            weight,
        })
    }

    fn is_valid(mask: &str) -> bool {
        let seed_re = Regex::new("^1+0[01]*1$").unwrap();
        seed_re.is_match(mask)
    }

    fn parse(mask: &str) -> Vec<u8> {
        mask.as_bytes()
            .iter()
            .flat_map(|b| match b {
                b'1' => Some(1),
                b'0' => Some(0),
                _ => None,
            })
            .collect()
    }

    fn get_differences(positions: &[usize]) -> Vec<usize> {
        // mask "1001101":
        //   positions    [0, 3, 4, 6]
        //   indices      [0, 1, 2, 3]
        //   differences  [0, 2, 2, 3]
        positions.iter().enumerate().map(|(i, &m)| m - i).collect()
    }

    fn get_positions(bytes: &[u8]) -> Vec<usize> {
        bytes
            .iter()
            .enumerate()
            .filter_map(|(i, &b)| (b == 1).then_some(i))
            .collect()
    }
}

impl Display for SeedMask {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.mask)
    }
}

#[cfg(test)]
mod tests {
    use super::{read_text_length, Header, SeedMask, SuffixSortType, OUTFILE_VERSION};
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
                num_suffixes: 17,
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

    #[test]
    fn test_seed_mask_new() -> Result<()> {
        let mask = SeedMask::new("110110101")?;
        let expected = SeedMask {
            mask: "110110101".to_string(),
            bytes: vec![1, 1, 0, 1, 1, 0, 1, 0, 1],
            positions: vec![0, 1, 3, 4, 6, 8],
            differences: vec![0, 0, 1, 1, 2, 3],
            weight: 6,
        };
        assert_eq!(mask, expected);

        assert!(SeedMask::new("").is_err());
        assert!(SeedMask::new("0").is_err());
        assert!(SeedMask::new("01").is_err());
        assert!(SeedMask::new("10").is_err());
        assert!(SeedMask::new("11").is_err());
        assert!(SeedMask::new("101").is_ok());
        Ok(())
    }

    #[test]
    fn test_seed_mask_from_bytes() -> Result<()> {
        let bytes: Vec<u8> = vec![1, 1, 0, 1, 1, 0, 1, 0, 1];
        let mask = SeedMask::from_bytes(&bytes)?;
        let expected = SeedMask {
            mask: "110110101".to_string(),
            bytes,
            positions: vec![0, 1, 3, 4, 6, 8],
            differences: vec![0, 0, 1, 1, 2, 3],
            weight: 6,
        };
        assert_eq!(mask, expected);
        Ok(())
    }

    #[test]
    fn test_create_seed_mask() -> Result<()> {
        // No 0 at beginning
        let res = SeedMask::new("0101");
        assert!(res.is_err());

        // No 0 at end
        let res = SeedMask::new("1010");
        assert!(res.is_err());

        // Nothing other than 0/1
        let res = SeedMask::new("1021");
        assert!(res.is_err());

        // Must have a 0
        let res = SeedMask::new("1111");
        assert!(res.is_err());

        // Valid
        let res = SeedMask::new("101");
        assert!(res.is_ok());

        let expected = SeedMask {
            mask: "101".to_string(),
            bytes: vec![1, 0, 1],
            positions: vec![0, 2],
            differences: vec![0, 1],
            weight: 2,
        };
        assert_eq!(res.unwrap(), expected);

        let res = SeedMask::new("11101101101000011");
        assert!(res.is_ok());

        let expected = SeedMask {
            mask: "11101101101000011".to_string(),
            bytes: vec![1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1],
            positions: vec![0, 1, 2, 4, 5, 7, 8, 10, 15, 16],
            differences: vec![0, 0, 0, 1, 1, 2, 2, 3, 7, 7],
            weight: 10,
        };
        assert_eq!(res.unwrap(), expected);

        Ok(())
    }

    #[test]
    fn test_display_seed_mask() -> Result<()> {
        for val in &["101", "11101101101000011"] {
            let mask = SeedMask::new(val).unwrap();
            assert_eq!(format!("{mask}"), val.to_string());
        }
        Ok(())
    }

    #[test]
    fn test_valid_seed_mask() -> Result<()> {
        let valid = ["101", "1001", "1101", "10101", "1110110110100001"];
        for pattern in valid {
            assert!(SeedMask::is_valid(pattern));
        }

        let invalid = [
            "", "abc", "1", "11", "111", "0", "00", "0111", "11100", "1a01",
        ];
        for pattern in invalid {
            assert!(!SeedMask::is_valid(pattern));
        }

        Ok(())
    }

    #[test]
    fn test_seed_mask_positions() -> Result<()> {
        assert_eq!(SeedMask::get_positions(&[1, 0, 1]), [0, 2]);
        assert_eq!(SeedMask::get_positions(&[1, 1, 0, 1, 1]), [0, 1, 3, 4]);
        assert_eq!(SeedMask::get_positions(&[1, 0, 1, 1, 0, 1]), [0, 2, 3, 5]);
        Ok(())
    }

    #[test]
    fn test_seed_mask_difference() -> Result<()> {
        assert_eq!(SeedMask::get_differences(&[]), []);

        // "11011" -> [0, 1, 3, 4]
        //           - 0  1  2  3
        //            ------------
        //            [0, 0, 1, 1]
        assert_eq!(SeedMask::get_differences(&[0, 1, 3, 4]), [0, 0, 1, 1]);

        // "100001" -> [0, 5]
        //            - 0  1
        //            --------
        //             [0, 4]
        assert_eq!(SeedMask::get_differences(&[0, 5]), [0, 4]);

        // "1001101" -> [0, 3, 4, 6]
        //             - 0  1  2  3
        //             -------------
        //              [0, 2, 2, 3]
        assert_eq!(SeedMask::get_differences(&[0, 3, 4, 6]), [0, 2, 2, 3]);
        Ok(())
    }
}
