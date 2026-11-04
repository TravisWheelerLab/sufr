//! Positioned file reads and writes, and on-disk access to the text/SA/LCP

use anyhow::{bail, Result};
use rayon::prelude::*;
use std::{
    cmp::min,
    fmt::{Debug, Display},
    fs::File,
    hash::Hash,
    io,
    marker::PhantomData,
    mem,
    ops::{Add, AddAssign, Range, Sub},
    slice,
};

mod sealed {
    pub trait Sealed {}
    impl Sealed for u8 {}
    impl Sealed for u32 {}
    impl Sealed for u64 {}
}

/// An element of a _.sufr_ array: `u8` for the text, and `u32` or `u64`
/// for the suffix and LCP arrays.
///
/// The trait is sealed because the byte casts in this module rely on
/// `u8`, `u32` and `u64` having no padding and no invalid bit patterns.
pub(crate) trait Int:
    Debug
    + AddAssign
    + Add<Output = Self>
    + Sub<Output = Self>
    + Copy
    + Default
    + Display
    + Ord
    + Hash
    + Send
    + Sync
    + sealed::Sealed
{
    /// The value as a `usize`.
    fn to_usize(&self) -> usize;

    /// The value of `val`, truncated to `Self`.
    fn from_usize(val: usize) -> Self;
}

impl Int for u8 {
    fn to_usize(&self) -> usize {
        *self as usize
    }

    fn from_usize(val: usize) -> Self {
        val as Self
    }
}

impl Int for u32 {
    fn to_usize(&self) -> usize {
        *self as usize
    }

    fn from_usize(val: usize) -> Self {
        val as Self
    }
}

impl Int for u64 {
    fn to_usize(&self) -> usize {
        *self as usize
    }

    fn from_usize(val: usize) -> Self {
        val as Self
    }
}

/// Read access to one array of a _.sufr_ file: the text, the suffix
/// array, or the LCP array.
#[derive(Debug)]
pub(crate) struct FileAccess<T: Int> {
    file: File,

    start_position: u64,

    end_position: u64,

    _marker: PhantomData<fn() -> T>,
}

impl<T: Int> FileAccess<T> {
    /// Open the array of `num_elements` values starting at byte `start`
    /// of `filename`.
    pub fn new(filename: &str, start: u64, num_elements: usize) -> Result<Self> {
        let file = File::open(filename)?;
        let size = num_elements * mem::size_of::<T>();
        Ok(FileAccess {
            file,
            start_position: start,
            end_position: start + size as u64,
            _marker: PhantomData,
        })
    }

    /// An iterator over the values of the array.
    pub fn iter(&self) -> FileAccessIter<'_, T> {
        FileAccessIter {
            file_access: self,
            buffer: vec![],
            buffer_pos: 0,
            current_position: self.start_position,
            exhausted: false,
        }
    }

    /// The value at `pos`, or `None` past the end of the array.
    ///
    /// The value at `pos`, `None` past the end of the array.
    pub fn get(&self, pos: usize) -> Result<Option<T>> {
        let seek = self.start_position + (pos * mem::size_of::<T>()) as u64;
        if seek < self.end_position {
            let mut val = T::from_usize(0);
            read_exact_at(
                &self.file,
                as_bytes_mut(std::slice::from_mut(&mut val)),
                seek,
            )?;
            Ok(Some(val))
        } else {
            Ok(None)
        }
    }

    /// The values in `range`, or an error when it runs past the end of
    /// the array.
    pub fn get_range(&self, range: Range<usize>) -> Result<Vec<T>> {
        assert!(range.start <= range.end);
        let start = self.start_position as usize + (range.start * mem::size_of::<T>());
        let end = self.start_position as usize + (range.end * mem::size_of::<T>());
        let valid = self.start_position as usize..self.end_position as usize + 1;
        if valid.contains(&start) && valid.contains(&end) {
            let num_vals = range.len();
            let mut buffer: Vec<T> = vec![T::from_usize(0); num_vals];
            read_exact_at(&self.file, as_bytes_mut(&mut buffer), start as u64)?;
            Ok(buffer)
        } else {
            bail!("Invalid range: {range:?}")
        }
    }
}

/// An iterator over the values of a [`FileAccess`].
#[derive(Debug)]
pub(crate) struct FileAccessIter<'a, T: Int> {
    file_access: &'a FileAccess<T>,
    buffer: Vec<T>,
    buffer_pos: usize,

    /// The file position of the first byte not yet in `buffer`.
    current_position: u64,

    exhausted: bool,
}

impl<'a, T: Int> FileAccessIter<'a, T> {
    /// The most bytes read from the file at once: 16 MiB.
    const BUFFER_SIZE: usize = 2usize.pow(24);
}

impl<T: Int> Iterator for FileAccessIter<'_, T> {
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        if self.exhausted {
            None
        } else {
            if self.buffer.is_empty() || self.buffer_pos == self.buffer.len() {
                if self.current_position >= self.file_access.end_position {
                    self.exhausted = true;
                    return None;
                }

                let bytes_wanted = min(
                    Self::BUFFER_SIZE,
                    (self.file_access.end_position - self.current_position) as usize,
                );
                let num_vals = bytes_wanted / mem::size_of::<T>();
                self.buffer.resize(num_vals, T::from_usize(0));
                read_exact_at(
                    &self.file_access.file,
                    as_bytes_mut(&mut self.buffer),
                    self.current_position,
                )
                .unwrap();

                self.current_position += (num_vals * mem::size_of::<T>()) as u64;
                self.buffer_pos = 0;
            }

            let val = self.buffer.get(self.buffer_pos).copied();

            self.buffer_pos += 1;
            val
        }
    }
}

/// The bytes of `buffer`, in native byte order.
pub(crate) fn as_bytes_mut<T: Int>(buffer: &mut [T]) -> &mut [u8] {
    // TODO: rework - endianness-dependent
    // SAFETY: T is u8, u32 or u64, which have no padding and
    //         accept every bit pattern, and the byte length is
    //         the length of `buffer`
    unsafe {
        slice::from_raw_parts_mut(
            buffer.as_mut_ptr() as *mut _,
            std::mem::size_of_val(buffer),
        )
    }
}

/// The bytes of `vec`, in native byte order.
pub(crate) fn as_bytes<T: Int>(vec: &[T]) -> &[u8] {
    // TODO: rework - endianness-dependent
    // SAFETY: T is u8, u32 or u64, which have no padding, and
    //         the byte length is the length of `vec`
    unsafe {
        slice::from_raw_parts(vec.as_ptr() as *const _, std::mem::size_of_val(vec))
    }
}

// the positioned reads and writes below take an absolute offset
// and leave the file cursor alone, so several threads may use
// one file at once

/// Write from `buf` at `offset`, returning the number of bytes written,
/// which may be fewer than requested.
#[cfg(unix)]
fn write_at(file: &File, buf: &[u8], offset: u64) -> io::Result<usize> {
    use std::os::unix::fs::FileExt;

    file.write_at(buf, offset)
}

#[cfg(windows)]
fn write_at(file: &File, buf: &[u8], offset: u64) -> io::Result<usize> {
    use std::os::windows::fs::FileExt;

    file.seek_write(buf, offset)
}

#[cfg(unix)]
fn read_exact_at(file: &File, buf: &mut [u8], offset: u64) -> io::Result<()> {
    use std::os::unix::fs::FileExt;

    file.read_exact_at(buf, offset)
}

// the loop of std's unix read_exact_at, over seek_read
#[cfg(windows)]
fn read_exact_at(file: &File, mut buf: &mut [u8], mut offset: u64) -> io::Result<()> {
    use std::os::windows::fs::FileExt;

    while !buf.is_empty() {
        match file.seek_read(buf, offset) {
            Ok(0) => break,
            Ok(n) => {
                buf = &mut buf[n..];
                offset += n as u64;
            }
            Err(ref e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    if !buf.is_empty() {
        Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            "failed to fill whole buffer",
        ))
    } else {
        Ok(())
    }
}

/// Write all of `buf` at `offset`, in parallel 64 MB chunks.
pub(crate) fn write_all_at(file: &File, buf: &[u8], offset: u64) -> Result<()> {
    let chunk_len = 1 << 26;
    buf.par_chunks(chunk_len)
        .enumerate()
        .try_for_each(|(i, mut chunk)| -> Result<()> {
            let mut pos = offset + (i * chunk_len) as u64;
            while !chunk.is_empty() {
                let n = write_at(file, chunk, pos)?;
                if n == 0 {
                    bail!("Failed to write to output file");
                }
                chunk = &chunk[n..];
                pos += n as u64;
            }
            Ok(())
        })
}

#[cfg(test)]
mod tests {
    use super::{as_bytes, as_bytes_mut, FileAccess};
    use crate::format::Header;
    use anyhow::Result;
    use pretty_assertions::assert_eq;

    #[test]
    fn test_as_bytes_mut() -> Result<()> {
        let mut buf: Vec<u32> = vec![1];
        as_bytes_mut(&mut buf).copy_from_slice(&[0, 0, 0, 0]);
        assert_eq!(buf, &[0u32]);

        let mut buf: Vec<u64> = vec![1];
        as_bytes_mut(&mut buf).copy_from_slice(&[0, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(buf, &[0u64]);

        let mut buf: Vec<u32> = vec![1];
        as_bytes_mut(&mut buf).copy_from_slice(&[1, 0, 0, 0]);
        assert_eq!(buf, &[1u32]);

        let mut buf: Vec<u64> = vec![1];
        as_bytes_mut(&mut buf).copy_from_slice(&[1, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(buf, &[1u64]);

        let mut buf: Vec<u32> = vec![1];
        as_bytes_mut(&mut buf).copy_from_slice(&[255, 255, 255, 255]);
        assert_eq!(buf, &[u32::MAX]);

        let mut buf: Vec<u64> = vec![1];
        as_bytes_mut(&mut buf)
            .copy_from_slice(&[255, 255, 255, 255, 255, 255, 255, 255]);
        assert_eq!(buf, &[u64::MAX]);

        let mut buf: Vec<u32> = vec![1, 1];
        as_bytes_mut(&mut buf).copy_from_slice(&[0, 0, 0, 0, 255, 255, 255, 255]);
        assert_eq!(buf, &[0u32, u32::MAX]);

        let mut buf: Vec<u64> = vec![1, 1];
        as_bytes_mut(&mut buf).copy_from_slice(&[
            0, 0, 0, 0, 0, 0, 0, 0, 255, 255, 255, 255, 255, 255, 255, 255,
        ]);
        assert_eq!(buf, &[0u64, u64::MAX]);

        Ok(())
    }

    #[test]
    fn test_as_bytes() -> Result<()> {
        let res = as_bytes(&[0u32]);
        assert_eq!(res, &[0, 0, 0, 0]);

        let res = as_bytes(&[0u64]);
        assert_eq!(res, &[0, 0, 0, 0, 0, 0, 0, 0]);

        let res = as_bytes(&[1u32]);
        assert_eq!(res, &[1, 0, 0, 0]);

        let res = as_bytes(&[1u64]);
        assert_eq!(res, &[1, 0, 0, 0, 0, 0, 0, 0]);

        let res = as_bytes(&[u32::MAX]);
        assert_eq!(res, &[255, 255, 255, 255]);

        let res = as_bytes(&[u64::MAX]);
        assert_eq!(res, &[255, 255, 255, 255, 255, 255, 255, 255]);

        let res = as_bytes(&[0u32, u32::MAX]);
        assert_eq!(res, &[0, 0, 0, 0, 255, 255, 255, 255]);

        let res = as_bytes(&[0u64, u64::MAX]);
        assert_eq!(
            res,
            &[0, 0, 0, 0, 0, 0, 0, 0, 255, 255, 255, 255, 255, 255, 255, 255]
        );

        Ok(())
    }

    #[test]
    fn test_file_access() -> Result<()> {
        let input_file = "../data/expected/abba.sufr";
        let header = Header::read(input_file, 4)?;
        let suffix_array_file: FileAccess<u32> = FileAccess::new(
            input_file,
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

        for (rank, &suffix) in suf_by_rank.iter().enumerate() {
            assert_eq!(suffix_array_file.get(rank)?, Some(suffix));
        }

        let res = suffix_array_file.get_range(1..100);
        assert!(res.is_err());
        assert_eq!(
            res.as_ref().unwrap_err().to_string(),
            "Invalid range: 1..100"
        );

        let res = suffix_array_file.get_range(8..9);
        assert!(res.is_ok());
        assert_eq!(res.as_ref().unwrap(), &[13]);

        let res = suffix_array_file.get_range(8..13);
        assert!(res.is_ok());
        assert_eq!(res.as_ref().unwrap(), &[13, 11, 9, 2, 4]);

        let res = suffix_array_file.get_range(1..8);
        assert!(res.is_ok());
        assert_eq!(res.as_ref().unwrap(), &[0, 12, 10, 1, 3, 5, 7]);

        let all: Vec<_> = suffix_array_file.iter().collect();
        assert_eq!(all, &[14, 0, 12, 10, 1, 3, 5, 7, 13, 11, 9, 2, 4, 6, 8]);

        for (i, suffix) in suffix_array_file.iter().enumerate() {
            assert_eq!(suf_by_rank[i], suffix);
        }

        Ok(())
    }
}
