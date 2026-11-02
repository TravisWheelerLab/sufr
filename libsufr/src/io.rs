//! Positioned file reads and writes, and on-disk access to the text/SA/LCP

use crate::types::Int;
use anyhow::{bail, Result};
use rayon::prelude::*;
use std::{cmp::min, fs::File, io, marker::PhantomData, mem, ops::Range, slice};

// --------------------------------------------------
/// Struct to mediate file access to on-disk arrays of text, suffix/LCP arrays
#[derive(Debug)]
pub struct FileAccess<T: Int> {
    /// A read-only filehandle to the _.sufr_ file
    file: File,

    /// The size in bytes for the entire text/SA/LCP
    pub size: usize,

    /// The starting byte position of the structure being read (text/SA/LCP)
    start_position: u64,

    /// The final byte position of the structure being read (text/SA/LCP)
    end_position: u64,

    _marker: PhantomData<fn() -> T>,
}

impl<T: Int> FileAccess<T> {
    /// Create a read-only file access to a portion of a _.sufr_ file
    /// representing the text, suffix array, or LCP array.
    /// This struct must be initialized using an `Int` of `u8` for the `text`
    /// or `u32`/`u64` for the SA/LCP.
    /// The metadata needed to create this can be found in the `Header`.
    ///
    /// Args:
    /// * `filename`: the _.sufr_ filename
    /// * `start`: the byte position in the file of the array
    /// * `num_elements`: the length of the text/SA/LCP
    pub fn new(filename: &str, start: u64, num_elements: usize) -> Result<Self> {
        let file = File::open(filename)?;
        let size = num_elements * mem::size_of::<T>();
        Ok(FileAccess {
            file,
            size,
            start_position: start,
            end_position: start + size as u64,
            _marker: PhantomData,
        })
    }

    /// Create a `FileAccessIter` iterator.
    pub fn iter(&self) -> FileAccessIter<'_, T> {
        FileAccessIter {
            file_access: self,
            buffer: vec![],
            buffer_pos: 0,
            current_position: self.start_position,
            exhausted: false,
        }
    }

    // --------------------------------------------------
    /// Return a value (`u8`/character from text or a SA/LCP value)
    ///
    /// Args:
    /// * `pos`: position in the array
    //
    // TODO: Ignoring lots of Results to return Option
    pub fn get(&self, pos: usize) -> Option<T> {
        // Don't bother looking for something beyond the end
        let seek = self.start_position + (pos * mem::size_of::<T>()) as u64;
        if seek < self.end_position {
            let mut val = T::from_usize(0);
            read_exact_at(
                &self.file,
                slice_int_to_slice_u8(std::slice::from_mut(&mut val)),
                seek,
            )
            .unwrap();
            Some(val)
        } else {
            None
        }
    }

    // --------------------------------------------------
    /// Return a range of values (`u8`/characters from text or a SA/LCP values)
    ///
    /// Args:
    /// * `range`: start/stop positions in the array
    pub fn get_range(&self, range: Range<usize>) -> Result<Vec<T>> {
        assert!(range.start <= range.end);
        let start = self.start_position as usize + (range.start * mem::size_of::<T>());
        let end = self.start_position as usize + (range.end * mem::size_of::<T>());
        let valid = self.start_position as usize..self.end_position as usize + 1;
        if valid.contains(&start) && valid.contains(&end) {
            let num_vals = range.len();
            let mut buffer: Vec<T> = vec![T::from_usize(0); num_vals];
            read_exact_at(
                &self.file,
                slice_int_to_slice_u8(&mut buffer),
                start as u64,
            )?;
            Ok(buffer)
        } else {
            bail!("Invalid range: {range:?}")
        }
    }
}

// --------------------------------------------------
/// An iterator over the values from a `FileAccess`
#[derive(Debug)]
pub struct FileAccessIter<'a, T: Int> {
    file_access: &'a FileAccess<T>,
    /// Internal buffer for reading a portion of the file
    buffer: Vec<T>,
    /// The current position when reading the buffer
    buffer_pos: usize,
    /// The current position after reading a portion of the structure
    /// from disk and into the `buffer`
    current_position: u64,
    /// Whether or not the last read off disk made it to the end of the array
    exhausted: bool,
}

impl<'a, T: Int> FileAccessIter<'a, T> {
    /// The maximum size in bytes of the buffer (currently 2^24, 16 MiB)
    const BUFFER_SIZE: usize = 2usize.pow(24);
}

impl<T: Int> Iterator for FileAccessIter<'_, T> {
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        if self.exhausted {
            None
        } else {
            // Fill the buffer
            if self.buffer.is_empty() || self.buffer_pos == self.buffer.len() {
                if self.current_position >= self.file_access.end_position {
                    self.exhausted = true;
                    return None;
                }

                // Read whole elements, at most BUFFER_SIZE bytes, straight
                // into the reused typed buffer
                let bytes_wanted = min(
                    Self::BUFFER_SIZE,
                    (self.file_access.end_position - self.current_position) as usize,
                );
                let num_vals = bytes_wanted / mem::size_of::<T>();
                self.buffer.resize(num_vals, T::from_usize(0));
                read_exact_at(
                    &self.file_access.file,
                    slice_int_to_slice_u8(&mut self.buffer),
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

// --------------------------------------------------
/// Convert a (mutable) slice of T into a
/// `&mut [u8]` (where `T` is the `Int` 32/64)
///
/// Args:
/// * `buffer`: slice of `T` values
pub(crate) fn slice_int_to_slice_u8<T: Int>(buffer: &mut [T]) -> &mut [u8] {
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

// --------------------------------------------------
/// Convert a `Vec<T>` (where `T` is the `Int` 32/64) into a
/// slice of raw U8 for serializing to disk
///
/// Args:
/// * `vec`: a vector of `T` values
pub(crate) fn vec_to_slice_u8<T: Int>(vec: &[T]) -> &[u8] {
    // TODO: rework - endianness-dependent
    // SAFETY: T is u8, u32 or u64, which have no padding, and
    //         the byte length is the length of `vec`
    unsafe {
        slice::from_raw_parts(
            vec.as_ptr() as *const _,
            std::mem::size_of_val(vec), //vec.len() * mem::size_of::<T>(),
        )
    }
}

// Position-independent reads and writes for both Unix and Windows

/// Write from `buf` at an absolute `offset` without using the file cursor,
/// so several threads may write disjoint regions of one file at once.
/// Returns the number of bytes written, which may be fewer than requested.
#[cfg(unix)]
fn write_at(file: &File, buf: &[u8], offset: u64) -> io::Result<usize> {
    use std::os::unix::fs::FileExt;

    file.write_at(buf, offset)
}

// Windows version of `write_at`; `seek_write` also takes an absolute
// offset per call, so concurrent calls place their data correctly.
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

// Essentially the same as the Unix read_exact_at implementation, but using seek_read
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

// --------------------------------------------------
#[cfg(test)]
mod tests {
    use super::{slice_int_to_slice_u8, vec_to_slice_u8, FileAccess};
    use crate::format::Header;
    use anyhow::Result;
    use pretty_assertions::assert_eq;

    #[test]
    fn test_slice_int_to_slice_u8() -> Result<()> {
        let mut buf: Vec<u32> = vec![1];
        slice_int_to_slice_u8(&mut buf).copy_from_slice(&[0, 0, 0, 0]);
        assert_eq!(buf, &[0u32]);

        let mut buf: Vec<u64> = vec![1];
        slice_int_to_slice_u8(&mut buf).copy_from_slice(&[0, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(buf, &[0u64]);

        let mut buf: Vec<u32> = vec![1];
        slice_int_to_slice_u8(&mut buf).copy_from_slice(&[1, 0, 0, 0]);
        assert_eq!(buf, &[1u32]);

        let mut buf: Vec<u64> = vec![1];
        slice_int_to_slice_u8(&mut buf).copy_from_slice(&[1, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(buf, &[1u64]);

        let mut buf: Vec<u32> = vec![1];
        slice_int_to_slice_u8(&mut buf).copy_from_slice(&[255, 255, 255, 255]);
        assert_eq!(buf, &[u32::MAX]);

        let mut buf: Vec<u64> = vec![1];
        slice_int_to_slice_u8(&mut buf)
            .copy_from_slice(&[255, 255, 255, 255, 255, 255, 255, 255]);
        assert_eq!(buf, &[u64::MAX]);

        let mut buf: Vec<u32> = vec![1, 1];
        slice_int_to_slice_u8(&mut buf)
            .copy_from_slice(&[0, 0, 0, 0, 255, 255, 255, 255]);
        assert_eq!(buf, &[0u32, u32::MAX]);

        let mut buf: Vec<u64> = vec![1, 1];
        slice_int_to_slice_u8(&mut buf).copy_from_slice(&[
            0, 0, 0, 0, 0, 0, 0, 0, 255, 255, 255, 255, 255, 255, 255, 255,
        ]);
        assert_eq!(buf, &[0u64, u64::MAX]);

        Ok(())
    }

    #[test]
    fn test_vec_to_slice_u8() -> Result<()> {
        let res = vec_to_slice_u8(&[0u32]);
        assert_eq!(res, &[0, 0, 0, 0]);

        let res = vec_to_slice_u8(&[0u64]);
        assert_eq!(res, &[0, 0, 0, 0, 0, 0, 0, 0]);

        let res = vec_to_slice_u8(&[1u32]);
        assert_eq!(res, &[1, 0, 0, 0]);

        let res = vec_to_slice_u8(&[1u64]);
        assert_eq!(res, &[1, 0, 0, 0, 0, 0, 0, 0]);

        let res = vec_to_slice_u8(&[u32::MAX]);
        assert_eq!(res, &[255, 255, 255, 255]);

        let res = vec_to_slice_u8(&[u64::MAX]);
        assert_eq!(res, &[255, 255, 255, 255, 255, 255, 255, 255]);

        let res = vec_to_slice_u8(&[0u32, u32::MAX]);
        assert_eq!(res, &[0, 0, 0, 0, 255, 255, 255, 255]);

        let res = vec_to_slice_u8(&[0u64, u64::MAX]);
        assert_eq!(
            res,
            &[0, 0, 0, 0, 0, 0, 0, 0, 255, 255, 255, 255, 255, 255, 255, 255]
        );

        Ok(())
    }

    // --------------------------------------------------
    #[test]
    fn test_file_access() -> Result<()> {
        let input_file = "../data/expected/abba.sufr";
        let header = Header::read(input_file, 4)?;
        let suffix_array_file: FileAccess<u32> = FileAccess::new(
            input_file,
            header.suffix_array_pos as u64,
            header.len_suffixes,
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
            let res = suffix_array_file.get(rank);
            assert!(res.is_some());
            assert_eq!(res.unwrap(), suffix);
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
