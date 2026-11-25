//! The fixtures every benchmark runs on, written once into `runs/`.
//!
//! Every byte is drawn from the generator. A pool of pre-drawn residues
//! recycled down the file would put a period in the text, and a suffix sort
//! compares the bytes, so a period of any length is what it would measure.
//!
//! The pathological fixtures are the normal one with one thing added, so a
//! row against the normal row measures that one thing.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use rand::{RngExt, SeedableRng, rngs::StdRng};

/// Fixed, so every run measures the same bytes.
const SEED: u64 = 1;

const LINE_WIDTH: usize = 60;

/// Records in a DNA fixture, each the same length.
pub(crate) const DNA_RECORDS: usize = 24;

const ACGT: &[u8] = b"ACGT";

// uniform random bases share prefixes of about log4(n)
// bases, which is not what a genome looks like: about a
// tenth of a mammalian genome is copies of a few hundred
// repeat families, each a few hundred bases long and 10-20%
// diverged from the family's consensus. the copies are what
// make a suffix sort's merges deep, so the fixture has them
//
//     elements     families, each ELEMENT_LEN random bases
//     one copy per GAP_MEAN random bases, on average
//     each base of a copy is replaced with DIVERGENCE_PERCENT
//     probability
const ELEMENT_LEN: usize = 300;
const ELEMENTS: usize = 50;
const GAP_MEAN: usize = 2_700;
const DIVERGENCE_PERCENT: u64 = 15;

/// Lengths of the runs of N in the long-n fixture, as a percentage of the
/// fixture, one run each.
//
// with every mitigation in libsufr off, a run costs about
// its length squared: at 10 MiB a run of 100,000 adds
// nothing, 300,000 adds 4x and 1,000,000 adds 34x
pub(crate) const N_RUN_PERCENTS: [usize; 3] = [1, 3, 10];

/// Copies of the unit in the long-repeats fixture.
pub(crate) const REPEAT_COPIES: usize = 10;

const AMINO: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

// log-normal record lengths, fitted to real proteomes:
//
//     median = exp(MU)             = 300 aa
//     mean   = exp(MU + SIGMA^2/2) ~ 383 aa
//     p99                          ~ 1530 aa
const MU: f64 = 5.704;
const SIGMA: f64 = 0.7;
const MIN_LEN: usize = 10;

/// Length of titin, clamping a freak draw from the tail to something a real
/// proteome could hold.
const MAX_LEN: usize = 35_000;

/// The normal DNA generator: uniform bases with diverged copies of a few
/// repeat families.
struct Background {
    rng: StdRng,
    elements: Vec<Vec<u8>>,
}

impl Background {
    fn new() -> Background {
        let mut rng = StdRng::seed_from_u64(SEED);
        let elements = (0..ELEMENTS)
            .map(|_| {
                (0..ELEMENT_LEN)
                    .map(|_| ACGT[rng.random_range(0..4)])
                    .collect()
            })
            .collect();

        Background { rng, elements }
    }

    /// Append `n` bases to `into`.
    fn fill(&mut self, into: &mut Vec<u8>, n: usize) {
        let want = into.len() + n;
        while into.len() < want {
            let gap = self.rng.random_range(0..=2 * GAP_MEAN);
            random_bases(&mut self.rng, into, gap);

            let element = &self.elements[self.rng.random_range(0..ELEMENTS)];
            for &base in element {
                into.push(if self.rng.random_range(0..100) < DIVERGENCE_PERCENT {
                    ACGT[self.rng.random_range(0..4)]
                } else {
                    base
                });
            }
        }
        into.truncate(want);
    }
}

/// Append `n` uniform bases to `into`, two bits of one draw per base.
fn random_bases(rng: &mut StdRng, into: &mut Vec<u8>, n: usize) {
    let mut left = n;
    while left > 0 {
        let word = rng.random::<u64>();
        for i in 0..32.min(left) {
            into.push(ACGT[((word >> (2 * i)) & 3) as usize]);
        }
        left -= 32.min(left);
    }
}

/// Write the normal DNA fixture of `bases` bases to `fasta`.
pub(crate) fn dna(fasta: &Path, bases: usize) -> std::io::Result<()> {
    let mut background = Background::new();
    let mut out = BufWriter::with_capacity(1 << 20, File::create(fasta)?);

    let len = bases / DNA_RECORDS;
    let mut record = Vec::with_capacity(len + 2 * GAP_MEAN + ELEMENT_LEN);

    for n in 0..DNA_RECORDS {
        record.clear();
        background.fill(&mut record, len);
        write_record(&mut out, &format!("chr{n}"), &record)?;
    }

    out.flush()
}

/// Write the long-n fixture to `fasta`: the normal DNA of `bases` bases
/// with one run of N inserted into the middle of each of the first records,
/// one run per entry of [`N_RUN_PERCENTS`].
pub(crate) fn long_n(fasta: &Path, bases: usize) -> std::io::Result<()> {
    let mut background = Background::new();
    let mut out = BufWriter::with_capacity(1 << 20, File::create(fasta)?);

    let len = bases / DNA_RECORDS;
    let mut record = Vec::with_capacity(len + 2 * GAP_MEAN + ELEMENT_LEN);

    for n in 0..DNA_RECORDS {
        record.clear();
        background.fill(&mut record, len);

        // inserted rather than written over, so the bases
        // around the run are the normal fixture's
        if let Some(percent) = N_RUN_PERCENTS.get(n) {
            let run = std::iter::repeat_n(b'N', bases * percent / 100);
            record.splice(len / 2..len / 2, run);
        }

        write_record(&mut out, &format!("chr{n}"), &record)?;
    }

    out.flush()
}

/// Write the long-repeats fixture to `fasta`: one record holding
/// [`REPEAT_COPIES`] copies of `bases / REPEAT_COPIES` normal bases, with
/// one random base between copies.
pub(crate) fn long_repeats(fasta: &Path, bases: usize) -> std::io::Result<()> {
    let mut background = Background::new();
    let mut out = BufWriter::with_capacity(1 << 20, File::create(fasta)?);

    let len = bases / REPEAT_COPIES;
    let mut unit = Vec::with_capacity(len + 2 * GAP_MEAN + ELEMENT_LEN);
    background.fill(&mut unit, len);

    let mut record = Vec::with_capacity(bases + REPEAT_COPIES);
    for n in 0..REPEAT_COPIES {
        if n > 0 {
            record.push(ACGT[background.rng.random_range(0..4)]);
        }
        record.extend_from_slice(&unit);
    }

    write_record(&mut out, "repeats", &record)?;
    out.flush()
}

/// Write the protein fixture of about `residues` residues to `fasta`.
pub(crate) fn protein(fasta: &Path, residues: usize) -> std::io::Result<()> {
    let mut rng = StdRng::seed_from_u64(SEED);
    let mut out = BufWriter::with_capacity(1 << 20, File::create(fasta)?);

    let mut written = 0;
    let mut records = 0;
    let mut record = Vec::with_capacity(MAX_LEN);

    while written < residues {
        let len = lognormal(&mut rng);
        record.clear();

        // a byte modulo 20 favours the first 16 residues by
        // 1/256 each, which a sort over them cannot see
        while record.len() < len {
            for byte in rng.random::<u64>().to_le_bytes() {
                record.push(AMINO[byte as usize % AMINO.len()]);
            }
        }
        record.truncate(len);

        write_record(&mut out, &format!("prot{records}"), &record)?;
        written += len;
        records += 1;
    }

    out.flush()
}

fn write_record(
    out: &mut impl Write,
    name: &str,
    record: &[u8],
) -> std::io::Result<()> {
    writeln!(out, ">{name}")?;
    for line in record.chunks(LINE_WIDTH) {
        out.write_all(line)?;
        out.write_all(b"\n")?;
    }
    Ok(())
}

/// One draw from the log-normal.
fn lognormal(rng: &mut StdRng) -> usize {
    // ln(0) is -inf, so a u1 of exactly zero has to be nudged
    let u1: f64 = rng.random::<f64>().max(f64::MIN_POSITIVE);
    let u2: f64 = rng.random::<f64>();

    // Box-Muller:
    //     z   = sqrt(-2 ln u1) * cos(2 pi u2)
    //     len = exp(MU + SIGMA * z)
    let z = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();

    ((MU + SIGMA * z).exp() as usize).clamp(MIN_LEN, MAX_LEN)
}
