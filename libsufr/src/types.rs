use std::{
    cmp::Ordering,
    fmt::{Debug, Display},
    hash::Hash,
    ops::Range,
    ops::{Add, Div, Sub},
};

// --------------------------------------------------
pub const OUTFILE_VERSION: u8 = 4;
pub const SENTINEL_CHARACTER: u8 = b'$';

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
