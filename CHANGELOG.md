# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added

- Start `CHANGELOG.md` for changes made after `0.7.12`
- (potentially breaking) `SufrBuilder` is parameterized with a new type
  bound `B: ScratchBuffer`, with implementations `DiskScratchBuffer`
  (default) and `MemoryScratchBuffer`. This makes it
  possible to use memory instead of disk for `SufrBuilder::new`.
- Implement vectorized LCP calculation (AVX2 on x86-64), significantly
  improving build time on some inputs
- Implement an LCP cache, significantly improving build time on some inputs with very long
  duplications or tandem repeats
- Add `SufrBuilderArgs::write_lcp` and `sufr create --write-lcp`. The LCP
  array is now optional and off by default. `search`, `count`, `locate`,
  and `extract` do not read it; `list --show-lcp` fails on a file without
  one. The header's LCP offset is 0 when there is no LCP array, and
  `SufrMetadata::has_lcp` and `summarize` report whether one is present.
- Subsample the suffix array for a query-time max query length without an
  LCP array by streaming the suffix array in blocks of four million
  elements and comparing each suffix with its predecessor in the text

### Changed

- (breaking) Simplify the `Int` trait definition and bounds, subsuming `FromUsize`
- (breaking) Make `Int` unable to implement from other crates, due to guarantees required by sufr's `unsafe` code
- (breaking) Change the partitioning strategy to be based on a prefix of each suffix,
  instead of using random sampling
- Reduce some unnecessary internal allocations
- During build, the merge sort makes better use of the available threads
  with work-stealing via `rayon::join`
- Write the output file with positioned writes (`pwrite` on Unix,
  `seek_write` on Windows) in parallel 64 MiB chunks instead of one
  sequential `BufWriter`, and patch the header with the array offsets at the end
- `FileAccessIter` refills its buffer with 16 MiB (2^24 bytes) of whole
  elements instead of 1 GiB, so iterating over an on-disk array smaller
  than 1 GiB no longer reads it whole

### Removed

- Remove `SufrBuilderArgs::random_seed` aka `sufr create -r,--random-seed`,
  which is not used for the prefix-based partitioning
- Remove and make private the `unsafe` byte-reading utilities
  `slice_u8_to_vec`, `usize_to_bytes`, `vec_to_slice_u8`.

### Fixed

- Disk I/O is updated to reduce intermediate allocations and remove unaligned memory reads/writes
- Searches (count, extract, or locate) with a `max_query_len`/`--max-query-len` could create cache
  files containing extraneous garbage memory; they no longer do so. The extra data was never read,
  so existing cache files should still be readable/writable across versions.
