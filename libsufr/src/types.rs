use anyhow::{bail, Result};
use regex::Regex;
use std::{
    cmp::Ordering,
    fmt::{self, Debug, Display},
    hash::Hash,
    ops::Range,
    ops::{Add, AddAssign, Div, Sub},
};

// --------------------------------------------------
pub const OUTFILE_VERSION: u8 = 6;
pub const SENTINEL_CHARACTER: u8 = b'$';

// --------------------------------------------------
#[derive(Debug)]
pub enum SuffixSortType {
    MaxQueryLen(usize), // can be zero
    Mask(SeedMask),
}

// --------------------------------------------------
#[derive(Debug, PartialEq)]
pub struct SeedMask {
    pub mask: String,
    pub bytes: Vec<u8>,
    pub positions: Vec<usize>,
    pub differences: Vec<usize>,
    pub weight: usize,
}

// --------------------------------------------------
impl SeedMask {
    pub fn new(mask: &str) -> Result<Self> {
        if !Self::valid(mask) {
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
        let positions = Self::get_positions(&bytes);
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

    pub fn valid(mask: &str) -> bool {
        let seed_re = Regex::new("^1+0[01]*1$").unwrap();
        seed_re.is_match(mask)
    }

    pub fn parse(mask: &str) -> Vec<u8> {
        mask.as_bytes()
            .iter()
            .flat_map(|b| match b {
                b'1' => Some(1),
                b'0' => Some(0),
                _ => None,
            })
            .collect()
    }

    pub fn get_differences(positions: &[usize]) -> Vec<usize> {
        // Mask: "1001101"
        // M: [0, 3, 4, 6]
        // U: [0, 1, 2, 3]
        // D: [0, 2, 2, 3]
        positions.iter().enumerate().map(|(i, &m)| m - i).collect()
    }

    pub fn get_positions(bytes: &[u8]) -> Vec<usize> {
        // [1, 0, 1] -> [0, 2]
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

// --------------------------------------------------
#[derive(Debug)]
pub struct SearchOptions {
    pub queries: Vec<String>,
    pub max_query_len: Option<usize>,
    pub low_memory: bool,
    pub find_suffixes: bool,
}

// --------------------------------------------------
#[derive(Debug)]
pub struct SearchResult<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    pub query_num: usize,
    pub query: String,
    pub locations: Option<SearchResultLocations<T>>,
}

// --------------------------------------------------
#[derive(Debug)]
pub struct SearchResultLocations<T> {
    pub ranks: Range<usize>,
    pub suffixes: Vec<T>,
}

// --------------------------------------------------
#[derive(Debug)]
pub struct Comparison {
    pub cmp: Ordering,
    pub lcp: usize,
}

// --------------------------------------------------
#[derive(Debug)]
pub struct SequenceFileData {
    pub seq: Vec<u8>,
    pub start_positions: Vec<usize>,
    pub headers: Vec<String>,
}

// --------------------------------------------------
pub trait Int:
    Debug
    + AddAssign
    + Add<Output = Self>
    + Sub<Output = Self>
    + Div<Output = Self>
    + Copy
    + Default
    + Display
    + Ord
    + Hash
    + serde::ser::Serialize
{
    fn to_usize(&self) -> usize;
}

impl Int for u32 {
    fn to_usize(&self) -> usize {
        *self as usize
    }
}

impl Int for u64 {
    fn to_usize(&self) -> usize {
        *self as usize
    }
}

pub trait FromUsize<T> {
    fn from_usize(val: usize) -> T;
}

impl FromUsize<u32> for u32 {
    fn from_usize(val: usize) -> u32 {
        val as u32
    }
}

impl FromUsize<u64> for u64 {
    fn from_usize(val: usize) -> u64 {
        val as u64
    }
}

// --------------------------------------------------
#[derive(Debug)]
pub struct ExtractOptions {
    pub queries: Vec<String>,
    pub max_query_len: Option<usize>,
    pub low_memory: bool,
    pub prefix_len: Option<usize>,
    pub suffix_len: Option<usize>,
}

// --------------------------------------------------
#[derive(Debug, PartialEq)]
pub struct ExtractResult {
    pub query_num: usize,
    pub query: String,
    pub sequences: Vec<ExtractSequence>,
}

// --------------------------------------------------
#[derive(Debug, PartialEq)]
pub struct ExtractSequence {
    pub suffix: usize,
    pub rank: usize,
    pub sequence_name: String,
    pub sequence_start: usize,
    pub sequence_range: Range<usize>,
    pub suffix_offset: usize,
}

// --------------------------------------------------
#[derive(Debug, PartialEq)]
pub struct LocateResult<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    pub query_num: usize,
    pub query: String,
    pub positions: Vec<LocatePosition<T>>,
}

// --------------------------------------------------
#[derive(Debug, PartialEq)]
pub struct LocatePosition<T>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + serde::ser::Serialize,
{
    pub suffix: T,
    pub rank: usize,
    pub sequence_name: String,
    pub sequence_position: T,
}

#[cfg(test)]
mod tests {
    use super::SeedMask;
    use anyhow::Result;
    use pretty_assertions::assert_eq;

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
            assert!(SeedMask::valid(pattern));
        }

        let invalid = [
            "", "abc", "1", "11", "111", "0", "00", "0111", "11100", "1a01",
        ];
        for pattern in invalid {
            assert!(!SeedMask::valid(pattern));
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
        // Empty is not a failure
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
