//! An array of values in memory or on disk.
//!
//! A [`Buffer`] is appended to during construction and read during a
//! search. Which one a function holds is a type parameter, so the
//! in-memory search compiles to slice reads and the on-disk one to
//! positioned reads, and a caller picks the trade when it opens the file.

use crate::io::{as_bytes, as_bytes_mut, read_exact_at, write_all_at, Int};
use anyhow::{bail, Result};
use std::{
    borrow::Cow,
    cmp::min,
    fmt::Debug,
    fs::{File, OpenOptions},
    marker::PhantomData,
    mem::size_of,
    ops::Range,
    path::Path,
    slice,
};
use tempfile::{NamedTempFile, TempPath};

/// Values addressable by index, appended to in order.
pub trait Buffer<W: Int>: Debug + Send + Sync + Sized {
    /// An empty buffer to append to.
    fn empty() -> Result<Self>;

    /// The `len` values at byte `start` of the file at `path`.
    fn open(path: &Path, start: u64, len: usize) -> Result<Self>;

    /// Append `vals`.
    fn extend_from_slice(&mut self, vals: &[W]) -> Result<()>;

    /// Write any appended values still in memory to the backing store.
    fn flush(&mut self) -> Result<()>;

    /// The number of values.
    fn len(&self) -> usize;

    /// Whether there are no values.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The value at `i`, or `None` past the end.
    fn get(&self, i: usize) -> Result<Option<W>>;

    /// The values in `range`, borrowed when they are in memory.
    //
    // &self, not &mut self: a search shares one buffer across
    // the rayon pool, so a Read or Seek bound would put a Mutex
    // behind every on-disk read
    fn range(&self, range: Range<usize>) -> Result<Cow<'_, [W]>>;
}

// ---

/// Values in memory.
//
// a newtype rather than an impl on Vec<W>: the trait's get and
// extend_from_slice would otherwise shadow the slice methods on
// every Vec in a module that imports Buffer
#[derive(Debug)]
pub struct Memory<W>(pub Vec<W>);

impl<W: Int> Buffer<W> for Memory<W> {
    fn empty() -> Result<Self> {
        Ok(Memory(vec![]))
    }

    fn open(path: &Path, start: u64, len: usize) -> Result<Self> {
        let mut vals = vec![W::default(); len];
        read_exact_at(&File::open(path)?, as_bytes_mut(&mut vals), start)?;
        Ok(Memory(vals))
    }

    fn extend_from_slice(&mut self, vals: &[W]) -> Result<()> {
        self.0.extend_from_slice(vals);
        Ok(())
    }

    fn flush(&mut self) -> Result<()> {
        Ok(())
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

/// Values in a file: a region of a _.sufr_ file, or a temporary file
/// that is deleted with the buffer.
#[derive(Debug)]
pub struct Disk<W> {
    handle: Handle,

    /// The byte offset of the first value.
    start: u64,

    /// The number of values in the file.
    written: usize,

    /// Appended values not yet in the file.
    pending: Vec<W>,

    _marker: PhantomData<fn() -> W>,
}

#[derive(Debug)]
enum Handle {
    Open(File),

    /// A temporary file, created at the first flush and opened for
    /// each access.
    //
    // a DNA build appends to up to 32768 buckets at once, most of
    // them empty: creating a file per bucket up front costs a
    // measurable share of a build, and holding them open is more
    // files than a process is usually allowed
    Temp(Option<TempPath>),
}

impl<W: Int> Disk<W> {
    /// The number of pending values that triggers a flush.
    const FLUSH_AT: usize = 4096;

    fn with_file<R>(&self, read: impl FnOnce(&File) -> Result<R>) -> Result<R> {
        match &self.handle {
            Handle::Open(file) => read(file),
            Handle::Temp(Some(path)) => read(&File::open(path)?),
            Handle::Temp(None) => bail!("no values have been written"),
        }
    }

    fn byte(&self, i: usize) -> u64 {
        self.start + (i * size_of::<W>()) as u64
    }
}

impl<W: Int> Buffer<W> for Disk<W> {
    fn empty() -> Result<Self> {
        Ok(Disk {
            handle: Handle::Temp(None),
            start: 0,
            written: 0,
            pending: vec![],
            _marker: PhantomData,
        })
    }

    fn open(path: &Path, start: u64, len: usize) -> Result<Self> {
        // a handle per region, not a dup of one: every thread of a
        // search reads through it, and preads on one file object
        // contend on its refcount
        Ok(Disk {
            handle: Handle::Open(File::open(path)?),
            start,
            written: len,
            pending: vec![],
            _marker: PhantomData,
        })
    }

    fn extend_from_slice(&mut self, vals: &[W]) -> Result<()> {
        self.pending.extend_from_slice(vals);
        if self.pending.len() >= Self::FLUSH_AT {
            self.flush()?;
        }
        Ok(())
    }

    fn flush(&mut self) -> Result<()> {
        if self.pending.is_empty() {
            return Ok(());
        }
        let at = self.byte(self.written);
        match &mut self.handle {
            Handle::Open(file) => write_all_at(file, as_bytes(&self.pending), at)?,
            Handle::Temp(Some(path)) => {
                let file = OpenOptions::new().append(true).open(path)?;
                write_all_at(&file, as_bytes(&self.pending), at)?;
            }
            Handle::Temp(path @ None) => {
                let (file, created) = NamedTempFile::new()?.into_parts();
                write_all_at(&file, as_bytes(&self.pending), at)?;
                *path = Some(created);
            }
        }
        self.written += self.pending.len();
        self.pending.clear();

        // the capacity, times the number of buckets, is
        // what the low-memory build would otherwise hold
        self.pending.shrink_to_fit();
        Ok(())
    }

    fn len(&self) -> usize {
        self.written + self.pending.len()
    }

    fn get(&self, i: usize) -> Result<Option<W>> {
        if i >= self.len() {
            return Ok(None);
        }
        if i >= self.written {
            return Ok(Some(self.pending[i - self.written]));
        }
        let mut val = W::default();
        self.with_file(|file| {
            read_exact_at(file, as_bytes_mut(slice::from_mut(&mut val)), self.byte(i))?;
            Ok(())
        })?;
        Ok(Some(val))
    }

    fn range(&self, range: Range<usize>) -> Result<Cow<'_, [W]>> {
        if range.start > range.end || range.end > self.len() {
            bail!("Invalid range: {range:?}");
        }
        let on_disk = range.start..min(range.end, self.written);
        let mut vals = vec![W::default(); on_disk.len()];
        if !on_disk.is_empty() {
            self.with_file(|file| {
                read_exact_at(file, as_bytes_mut(&mut vals), self.byte(on_disk.start))?;
                Ok(())
            })?;
        }
        if range.end > self.written {
            let from = range.start.saturating_sub(self.written);
            vals.extend_from_slice(&self.pending[from..range.end - self.written]);
        }
        Ok(Cow::Owned(vals))
    }
}

#[cfg(test)]
mod tests {
    use super::{Buffer, Disk, Memory};
    use crate::format::Header;
    use anyhow::Result;
    use pretty_assertions::assert_eq;
    use std::path::Path;

    #[test]
    fn test_open_region() -> Result<()> {
        let input_file = "../data/expected/abba.sufr";
        let header = Header::read(input_file, 4)?;
        let file = Path::new(input_file);
        let suffixes: Disk<u32> =
            Disk::open(file, header.suffix_array_pos as u64, header.num_suffixes)?;
        let in_memory: Memory<u32> =
            Memory::open(file, header.suffix_array_pos as u64, header.num_suffixes)?;
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
        Ok(())
    }

    // the tail of a disk buffer is in memory until it is flushed, and a
    // read covers the file, the tail and a range over both
    #[test]
    fn test_append_and_read_back() -> Result<()> {
        let flush_at = Disk::<u32>::FLUSH_AT;
        let vals: Vec<u32> = (0..flush_at as u32 + 10).collect();
        let mut disk: Disk<u32> = Disk::empty()?;
        let mut mem: Memory<u32> = Memory::empty()?;
        for chunk in vals.chunks(flush_at / 2 + 1) {
            disk.extend_from_slice(chunk)?;
            mem.extend_from_slice(chunk)?;
        }
        assert_eq!(disk.len(), vals.len());
        assert_eq!(mem.len(), vals.len());
        assert_eq!(disk.get(0)?, Some(0));
        assert_eq!(disk.get(vals.len() - 1)?, Some(vals.len() as u32 - 1));
        assert_eq!(disk.get(vals.len())?, None);
        assert_eq!(&*disk.range(0..vals.len())?, &vals);
        assert_eq!(
            &*disk.range(flush_at - 2..flush_at + 2)?,
            &vals[flush_at - 2..flush_at + 2]
        );
        assert_eq!(
            &*disk.range(flush_at + 1..flush_at + 3)?,
            &vals[flush_at + 1..flush_at + 3]
        );
        assert_eq!(&*disk.range(5..5)?, &[]);

        disk.flush()?;
        assert_eq!(disk.len(), vals.len());
        assert_eq!(&*disk.range(0..vals.len())?, &vals);
        assert_eq!(&*mem.range(0..vals.len())?, &vals);
        Ok(())
    }

    #[test]
    fn test_temp_file_is_removed() -> Result<()> {
        let mut disk: Disk<u32> = Disk::empty()?;
        disk.extend_from_slice(&[1, 2, 3])?;
        disk.flush()?;
        let path = match &disk.handle {
            super::Handle::Temp(Some(path)) => path.to_path_buf(),
            _ => unreachable!(),
        };
        assert!(path.is_file());
        drop(disk);
        assert!(!path.exists());
        Ok(())
    }
}
