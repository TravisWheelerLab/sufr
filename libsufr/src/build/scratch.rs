use crate::io::{as_bytes, as_bytes_mut, Int};
use anyhow::Result;
use std::{
    borrow::Cow,
    fs::{self, File, OpenOptions},
    io::{Read, Write},
    path::PathBuf,
};
use tempfile::NamedTempFile;

/// Where a [`Scratch`] buffer keeps its values.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum ScratchKind {
    Disk,
    Memory,
}

/// A buffer of suffixes or LCPs during construction, in memory or in a
/// temporary file.
#[derive(Debug)]
pub(super) enum Scratch<T: Int> {
    Disk(DiskScratch<T>),
    Memory(Vec<T>),
}

impl<T: Int> Scratch<T> {
    /// An empty buffer of `kind`.
    pub(super) fn new(kind: ScratchKind) -> Self {
        match kind {
            ScratchKind::Disk => Scratch::Disk(DiskScratch {
                buf: vec![],
                count: 0,
                path: None,
            }),
            ScratchKind::Memory => Scratch::Memory(vec![]),
        }
    }

    pub(super) fn kind(&self) -> ScratchKind {
        match self {
            Scratch::Disk(_) => ScratchKind::Disk,
            Scratch::Memory(_) => ScratchKind::Memory,
        }
    }

    /// Append `vals`.
    pub(super) fn extend_from_slice(&mut self, vals: &[T]) -> Result<()> {
        match self {
            Scratch::Disk(disk) => disk.extend_from_slice(vals),
            Scratch::Memory(data) => {
                data.extend_from_slice(vals);
                Ok(())
            }
        }
    }

    /// Write any values still in memory to the backing store.
    pub(super) fn flush(&mut self) -> Result<()> {
        match self {
            Scratch::Disk(disk) => disk.flush(),
            Scratch::Memory(_) => Ok(()),
        }
    }

    /// The number of values appended.
    pub(super) fn count(&self) -> usize {
        match self {
            Scratch::Disk(disk) => disk.count + disk.buf.len(),
            Scratch::Memory(data) => data.len(),
        }
    }

    /// All the values, freeing the backing store.
    pub(super) fn consume(self) -> Result<Vec<T>> {
        match self {
            Scratch::Disk(disk) => disk.consume(),
            Scratch::Memory(data) => Ok(data),
        }
    }
}

/// A scratch buffer in a temporary file.
#[derive(Debug)]
pub(super) struct DiskScratch<T: Int> {
    buf: Vec<T>,
    count: usize,
    path: Option<PathBuf>,
}

impl<T: Int> DiskScratch<T> {
    const FLUSH_AT: usize = 4096;

    fn open(&mut self) -> Result<File> {
        if self.path.is_none() {
            let tmp = NamedTempFile::new()?;
            let (_, path) = tmp.keep()?;
            self.path = Some(path)
        }
        Ok(OpenOptions::new()
            .create(true)
            .append(true)
            .open(self.path.as_ref().unwrap())?)
    }

    fn extend_from_slice(&mut self, vals: &[T]) -> Result<()> {
        self.buf.extend_from_slice(vals);
        if self.buf.len() >= Self::FLUSH_AT {
            self.flush()?;
        }
        Ok(())
    }

    fn flush(&mut self) -> Result<()> {
        if !self.buf.is_empty() {
            let mut file = self.open()?;
            file.write_all(as_bytes(&self.buf))?;
            self.count += self.buf.len();
            self.buf.clear();
            self.buf.shrink_to_fit();
        }
        Ok(())
    }

    fn read(&self) -> Result<Cow<'_, [T]>> {
        match &self.path {
            Some(path) => {
                let mut file = File::open(path)?;
                let size = file.metadata()?.len() as usize;

                let file_count = size / std::mem::size_of::<T>();
                let mut data = Vec::with_capacity(file_count + self.buf.len());
                data.resize(file_count, T::from_usize(0));
                file.read_exact(as_bytes_mut(&mut data))?;

                data.extend_from_slice(&self.buf);
                Ok(Cow::Owned(data))
            }
            None => Ok(Cow::Borrowed(&self.buf)),
        }
    }

    fn consume(mut self) -> Result<Vec<T>> {
        let buf = match self.read()? {
            Cow::Owned(b) => b,
            Cow::Borrowed(b) => b.to_vec(),
        };
        if let Some(path) = self.path.take() {
            fs::remove_file(path)?;
        }
        Ok(buf)
    }
}

impl<T: Int> Drop for DiskScratch<T> {
    fn drop(&mut self) {
        if let Some(path) = self.path.take() {
            // best effort, since this is Drop and consume()
            // may have removed the file already
            let _ = fs::remove_file(path);
        }
    }
}
