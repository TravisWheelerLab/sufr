//! An array of values the search reads: a region of a _.sufr_ file, or
//! the same values in memory.
//!
//! Which one [`crate::typed::TypedSuffixArray`] holds for its text and
//! for its suffix array is a type parameter, so the in-memory search
//! compiles to slice reads and the on-disk one to positioned reads, and
//! a caller picks the trade when it opens the file.

use crate::io::{as_bytes_mut, read_exact_at, Int};
use anyhow::{bail, Result};
use std::{
    borrow::Cow, fmt::Debug, fs::File, marker::PhantomData, mem::size_of, ops::Range,
    path::Path, slice,
};

/// Values addressable by index.
pub trait Source<W: Int>: Debug + Send + Sync + Sized {
    /// The `len` values at byte `start` of the file at `path`.
    fn open(path: &Path, start: u64, len: usize) -> Result<Self>;

    /// The number of values.
    fn len(&self) -> usize;

    /// The value at `i`, or `None` past the end.
    fn get(&self, i: usize) -> Result<Option<W>>;

    /// The values in `range`, borrowed when they are in memory.
    //
    // &self, not &mut self: a search shares one source across
    // the rayon pool, so a Read or Seek bound would put a Mutex
    // behind every on-disk read
    fn range(&self, range: Range<usize>) -> Result<Cow<'_, [W]>>;
}

// ---

/// Values in memory.
#[derive(Debug)]
pub struct MemorySource<W>(pub Vec<W>);

impl<W: Int> Source<W> for MemorySource<W> {
    fn open(path: &Path, start: u64, len: usize) -> Result<Self> {
        let mut vals = vec![W::default(); len];
        read_exact_at(&File::open(path)?, as_bytes_mut(&mut vals), start)?;
        Ok(MemorySource(vals))
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn get(&self, i: usize) -> Result<Option<W>> {
        Ok(self.0.get(i).copied())
    }

    fn range(&self, range: Range<usize>) -> Result<Cow<'_, [W]>> {
        match self.0.get(range.clone()) {
            Some(vals) => Ok(Cow::Borrowed(vals)),
            None => bail!("Invalid range: {range:?}"),
        }
    }
}

// ---

/// Values in a region of a file, read as they are asked for.
#[derive(Debug)]
pub struct DiskSource<W> {
    /// A handle of its own, not a dup of one shared with the other
    /// regions.
    //
    // every thread of a search reads through it, and preads on
    // one file object contend on its refcount
    file: File,

    /// The byte offset of the first value.
    start: u64,

    len: usize,

    _marker: PhantomData<fn() -> W>,
}

impl<W: Int> DiskSource<W> {
    fn byte(&self, i: usize) -> u64 {
        self.start + (i * size_of::<W>()) as u64
    }
}

impl<W: Int> Source<W> for DiskSource<W> {
    fn open(path: &Path, start: u64, len: usize) -> Result<Self> {
        Ok(DiskSource {
            file: File::open(path)?,
            start,
            len,
            _marker: PhantomData,
        })
    }

    fn len(&self) -> usize {
        self.len
    }

    fn get(&self, i: usize) -> Result<Option<W>> {
        if i >= self.len {
            return Ok(None);
        }
        let mut val = W::default();
        read_exact_at(
            &self.file,
            as_bytes_mut(slice::from_mut(&mut val)),
            self.byte(i),
        )?;
        Ok(Some(val))
    }

    fn range(&self, range: Range<usize>) -> Result<Cow<'_, [W]>> {
        if range.start > range.end || range.end > self.len {
            bail!("Invalid range: {range:?}");
        }
        let mut vals = vec![W::default(); range.len()];
        read_exact_at(&self.file, as_bytes_mut(&mut vals), self.byte(range.start))?;
        Ok(Cow::Owned(vals))
    }
}

#[cfg(test)]
mod tests {
    use super::{DiskSource, MemorySource, Source};
    use crate::format::Header;
    use anyhow::Result;
    use pretty_assertions::assert_eq;
    use std::path::Path;

    #[test]
    fn test_open_region() -> Result<()> {
        let input_file = "../data/expected/abba.sufr";
        let header = Header::read(input_file, 4)?;
        let file = Path::new(input_file);
        let suffixes: DiskSource<u32> = DiskSource::open(
            file,
            header.suffix_array_pos as u64,
            header.num_suffixes,
        )?;
        let in_memory: MemorySource<u32> = MemorySource::open(
            file,
            header.suffix_array_pos as u64,
            header.num_suffixes,
        )?;
        let suf_by_rank = [
            14, //  0: #
            0,  //  1: AABABABABBABAB#
            12, //  2: AB#
            10, //  3: ABAB#
            1,  //  4: ABABABABBABAB#
            3,  //  5: ABABABBABAB#
            5,  //  6: ABABBABAB#
            7,  //  7: ABBABAB#
            13, //  8: B#
            11, //  9: BAB#
            9,  // 10: BABAB#
            2,  // 11: BABABABBABAB#
            4,  // 12: BABABBABAB#
            6,  // 13: BABBABAB#
            8,  // 14: BBABAB#
        ];
        assert_eq!(in_memory.0, suf_by_rank);
        assert_eq!(suffixes.len(), 15);

        for (rank, &suffix) in suf_by_rank.iter().enumerate() {
            assert_eq!(suffixes.get(rank)?, Some(suffix));
            assert_eq!(in_memory.get(rank)?, Some(suffix));
        }
        assert_eq!(suffixes.get(15)?, None);
        assert_eq!(in_memory.get(15)?, None);

        let res = suffixes.range(1..100);
        assert_eq!(res.unwrap_err().to_string(), "Invalid range: 1..100");
        let res = in_memory.range(1..100);
        assert_eq!(res.unwrap_err().to_string(), "Invalid range: 1..100");

        assert_eq!(&*suffixes.range(8..9)?, &[13]);
        assert_eq!(&*suffixes.range(8..13)?, &[13, 11, 9, 2, 4]);
        assert_eq!(&*suffixes.range(1..8)?, &[0, 12, 10, 1, 3, 5, 7]);
        assert_eq!(&*suffixes.range(0..15)?, &suf_by_rank);
        assert_eq!(&*in_memory.range(1..8)?, &[0, 12, 10, 1, 3, 5, 7]);
        assert_eq!(&*in_memory.range(5..5)?, &[]);
        Ok(())
    }
}
