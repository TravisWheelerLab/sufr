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

### Changed

- (breaking) Simplify the `Int` trait definition and bounds, subsuming `FromUsize`
- (breaking) Change the partitioning strategy to be based on a prefix of each suffix,
  instead of using random sampling
- Reduce some unnecessary internal allocations
- During build, the merge sort makes better use of the available threads
  with work-stealing via `rayon::join`

### Removed

- Remove `SufrBuilderArgs::random_seed` aka `sufr create -r,--random-seed`,
  which is not used for the prefix-based partitioning
