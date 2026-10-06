//! Tests that subsampling a suffix array for a shorter max query length
//! gives the same result with and without an LCP array on disk.

use crate::{sufr_builder::SufrBuilder, sufr_file::SufrFile, types::SufrBuilderArgs};
use anyhow::Result;
use pretty_assertions::assert_eq;
use tempfile::NamedTempFile;

/// A small deterministic generator (xorshift64*) so the test needs no
/// random-number crate
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545F4914F6CDD1D)
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

fn random_text(rng: &mut Rng) -> Vec<u8> {
    let mut text = vec![];
    for rec in 0..3 {
        if rec > 0 {
            text.push(b'%');
        }
        let len = 300 + rng.below(400);
        text.extend((0..len).map(|_| b"ACGT"[rng.below(4)]));
        // A few Ns so some suffixes are skipped in DNA mode
        let at = rng.below(len - 5);
        let start = text.len() - len + at;
        for b in text[start..start + 4].iter_mut() {
            *b = b'N';
        }
    }
    text.push(b'$');
    text
}

/// Build with or without an LCP array and return the open file
fn build(
    text: &[u8],
    seed_mask: Option<&str>,
    max_query_len: Option<usize>,
    write_lcp: bool,
) -> Result<(NamedTempFile, SufrFile<u32>)> {
    let outfile = NamedTempFile::new()?;
    let outpath = outfile.path().to_string_lossy().to_string();
    SufrBuilder::<u32>::new(SufrBuilderArgs {
        text: text.to_vec(),
        path: Some(outpath.clone()),
        low_memory: false,
        max_query_len,
        is_dna: true,
        allow_ambiguity: false,
        ignore_softmask: false,
        sequence_starts: vec![0],
        sequence_names: vec!["s".to_string()],
        num_partitions: 3,
        seed_mask: seed_mask.map(String::from),
        write_lcp,
    })?;
    let sufr_file: SufrFile<u32> = SufrFile::read(&outpath, false)?;
    Ok((outfile, sufr_file))
}

/// Independent definition: keep rank i when the first `mql` key symbols of
/// SA[i] differ from those of SA[i-1] (a symbol past the text end differs
/// from everything, including another past-the-end symbol).
fn expected_subsample(
    text: &[u8],
    sa: &[u32],
    offsets: &[usize],
    mql: usize,
) -> (Vec<u32>, Vec<u32>) {
    let key = |pos: usize| -> Vec<Option<u8>> {
        offsets[..mql.min(offsets.len())]
            .iter()
            .map(|&off| text.get(pos + off).copied())
            .collect()
    };
    let mut sub_sa = vec![];
    let mut rank = vec![];
    for i in 0..sa.len() {
        let keep = i == 0 || {
            let a = key(sa[i - 1] as usize);
            let b = key(sa[i] as usize);
            !a.iter().zip(&b).all(|(x, y)| x.is_some() && x == y)
        };
        if keep {
            sub_sa.push(sa[i]);
            rank.push(i as u32);
        }
    }
    (sub_sa, rank)
}

#[test]
fn subsample_without_lcp_matches_definition_and_lcp_version() -> Result<()> {
    let mut rng = Rng(11);
    let text = random_text(&mut rng);

    // Seed-mask index, query-time MQL shorter than the weight
    let mask = "11101101101111";
    let offsets: Vec<usize> = mask
        .bytes()
        .enumerate()
        .filter(|(_, b)| *b == b'1')
        .map(|(i, _)| i)
        .collect();
    for mql in [1, 2, 4, 7, 10] {
        let (_f1, mut with_lcp) = build(&text, Some(mask), None, true)?;
        let (_f2, mut without_lcp) = build(&text, Some(mask), None, false)?;
        assert!(with_lcp.has_lcp);
        assert!(!without_lcp.has_lcp);
        let sa: Vec<u32> = without_lcp.suffix_array_file.iter().collect();
        let expected = expected_subsample(&text, &sa, &offsets, mql);
        let from_text = without_lcp.subsample_suffix_array(mql);
        assert_eq!(from_text, expected, "mask, mql={mql}, vs definition");
        let from_lcp = with_lcp.subsample_suffix_array(mql);
        assert_eq!(from_text, from_lcp, "mask, mql={mql}, vs LCP array");
    }

    // Max-query-len index, query-time MQL shorter than the build-time one
    let built_mql = 12;
    let offsets: Vec<usize> = (0..built_mql).collect();
    for mql in [1, 3, 6, 11] {
        let (_f1, mut with_lcp) = build(&text, None, Some(built_mql), true)?;
        let (_f2, mut without_lcp) = build(&text, None, Some(built_mql), false)?;
        let sa: Vec<u32> = without_lcp.suffix_array_file.iter().collect();
        let expected = expected_subsample(&text, &sa, &offsets, mql);
        let from_text = without_lcp.subsample_suffix_array(mql);
        assert_eq!(from_text, expected, "mql index, mql={mql}, vs definition");
        let from_lcp = with_lcp.subsample_suffix_array(mql);
        assert_eq!(from_text, from_lcp, "mql index, mql={mql}, vs LCP array");
    }
    Ok(())
}
