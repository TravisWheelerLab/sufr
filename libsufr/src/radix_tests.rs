//! Tests that the radix-sort construction paths agree with the merge-sort
//! path and with an independent definition of the seed-mask order.

use crate::{
    subsample_tests::Rng,
    sufr_builder::SufrBuilder,
    sufr_file::SufrFile,
    types::{CountOptions, SeedMask, SortStrategy, SufrBuilderArgs},
    util::read_sequence_file,
};
use anyhow::Result;
use pretty_assertions::assert_eq;
use std::path::Path;
use tempfile::NamedTempFile;

const STRATEGIES: [SortStrategy; 2] =
    [SortStrategy::RadixInMemory, SortStrategy::RadixPartitioned];

/// Inputs to one build
#[derive(Clone)]
struct Case {
    text: Vec<u8>,
    starts: Vec<usize>,
    names: Vec<String>,
    is_dna: bool,
    allow_ambiguity: bool,
    seed_mask: Option<String>,
    max_query_len: Option<usize>,
}

impl Case {
    fn args(
        &self,
        path: String,
        num_partitions: usize,
        sort_strategy: SortStrategy,
    ) -> SufrBuilderArgs {
        SufrBuilderArgs {
            text: self.text.clone(),
            path: Some(path),
            low_memory: true,
            max_query_len: self.max_query_len,
            is_dna: self.is_dna,
            allow_ambiguity: self.allow_ambiguity,
            ignore_softmask: false,
            sequence_starts: self.starts.clone(),
            sequence_names: self.names.clone(),
            num_partitions,
            seed_mask: self.seed_mask.clone(),
            sort_strategy,
            write_lcp: true,
        }
    }

    fn describe(&self) -> String {
        format!(
            "len={} dna={} ambig={} mask={:?} mql={:?}",
            self.text.len(),
            self.is_dna,
            self.allow_ambiguity,
            self.seed_mask,
            self.max_query_len
        )
    }
}

/// Build with a strategy and return (SA, LCP)
fn build(
    case: &Case,
    strategy: SortStrategy,
    num_partitions: usize,
) -> Result<(Vec<u32>, Vec<u32>)> {
    let outfile = NamedTempFile::new()?;
    let outpath = outfile.path().to_string_lossy().to_string();
    SufrBuilder::<u32>::new(case.args(outpath.clone(), num_partitions, strategy))?;
    let sufr_file: SufrFile<u32> = SufrFile::read(&outpath, false)?;
    let sa: Vec<u32> = sufr_file.suffix_array_file.iter().collect();
    let lcp: Vec<u32> = sufr_file.lcp_file.iter().collect();
    Ok((sa, lcp))
}

fn assert_strategies_agree(case: &Case, num_partitions: usize) -> Result<()> {
    let (sa, lcp) = build(case, SortStrategy::Merge, 1)?;
    for strategy in STRATEGIES {
        let (sa2, lcp2) = build(case, strategy, num_partitions)?;
        assert_eq!(sa2, sa, "SA differs for {strategy:?} on {}", case.describe());
        assert_eq!(lcp2, lcp, "LCP differs for {strategy:?} on {}", case.describe());
    }
    Ok(())
}

/// The key offsets of a case: the mask's care positions or `0..mql`
fn key_offsets(case: &Case) -> Vec<usize> {
    match (&case.seed_mask, case.max_query_len) {
        (Some(mask), _) => SeedMask::new(mask).unwrap().positions,
        (None, Some(mql)) => (0..mql).collect(),
        _ => panic!("case needs a mask or an MQL"),
    }
}

/// Bytes at the key offsets; `None` past the end of the text
fn masked(text: &[u8], offsets: &[usize], pos: usize) -> Vec<Option<u8>> {
    offsets.iter().map(|&off| text.get(pos + off).copied()).collect()
}

/// Independent definition of the expected order: indexed positions sorted
/// by the bytes at the key offsets, with "past the end" sorting first and
/// ties broken by descending position.
fn expected_order(case: &Case) -> Vec<usize> {
    let text = &case.text;
    let offsets = key_offsets(case);
    let mut expected: Vec<usize> = (0..text.len())
        .filter(|&i| {
            text[i] == b'$'
                || !case.is_dna
                || b"ACGT".contains(&text[i])
                || case.allow_ambiguity
        })
        .collect();
    expected.sort_by(|&a, &b| {
        masked(text, &offsets, a)
            .cmp(&masked(text, &offsets, b))
            .then(b.cmp(&a))
    });
    expected
}

/// Expected LCP (in key symbols) of two positions
fn expected_lcp(case: &Case, a: usize, b: usize) -> usize {
    let offsets = key_offsets(case);
    masked(&case.text, &offsets, a)
        .iter()
        .zip(&masked(&case.text, &offsets, b))
        .take_while(|(x, y)| x.is_some() && x == y)
        .count()
}

fn assert_matches_definition(case: &Case, num_partitions: usize) -> Result<()> {
    let expected = expected_order(case);
    for strategy in STRATEGIES {
        let (sa, lcp) = build(case, strategy, num_partitions)?;
        let sa: Vec<usize> = sa.into_iter().map(|v| v as usize).collect();
        assert_eq!(sa, expected, "{strategy:?} {}", case.describe());
        assert_eq!(lcp[0], 0);
        for k in 1..sa.len() {
            assert_eq!(
                lcp[k] as usize,
                expected_lcp(case, sa[k - 1], sa[k]),
                "LCP at rank {k} for {strategy:?} {}",
                case.describe()
            );
        }
    }
    Ok(())
}

/// A value in `lo..hi`
fn between(rng: &mut Rng, lo: usize, hi: usize) -> usize {
    lo + rng.below(hi - lo)
}

fn random_dna(rng: &mut Rng, len: usize) -> Vec<u8> {
    (0..len).map(|_| b"ACGT"[rng.below(4)]).collect()
}

/// Several records with soft-masked bases, short N runs, and one long N
/// run; records are joined with `%` and the text ends in `$`.
fn random_genome(rng: &mut Rng, long_n_run: bool) -> (Vec<u8>, Vec<usize>, Vec<String>) {
    let mut text = vec![];
    let mut starts = vec![];
    let mut names = vec![];
    for rec in 0..3 {
        if rec > 0 {
            text.push(b'%');
        }
        starts.push(text.len());
        names.push(format!("seq{rec}"));
        let len = between(rng, 400, 900);
        let mut seq = random_dna(rng, len);
        for _ in 0..3 {
            let at = rng.below(len - 10);
            for b in seq[at..at + between(rng, 1, 10)].iter_mut() {
                *b = b'N';
            }
        }
        if rec == 1 && long_n_run {
            seq.extend(vec![b'N'; 1200]);
            seq.extend(random_dna(rng, 300));
        }
        text.extend(seq);
    }
    text.push(b'$');
    (text, starts, names)
}

fn mask_cases(text: Vec<u8>, starts: Vec<usize>, names: Vec<String>) -> Vec<Case> {
    let mut cases = vec![];
    for mask in ["101", "11011", "11000111", "1101101", "11101101101111"] {
        for (is_dna, allow_ambiguity) in [(true, false), (true, true), (false, false)] {
            cases.push(Case {
                text: text.clone(),
                starts: starts.clone(),
                names: names.clone(),
                is_dna,
                allow_ambiguity,
                seed_mask: Some(mask.to_string()),
                max_query_len: None,
            });
        }
    }
    cases
}

// --------------------------------------------------
#[test]
fn radix_agrees_with_merge_on_test_inputs() -> Result<()> {
    for file in ["1.fa", "2.fa", "3.fa", "mostlya1.fa", "mostlya2.fa", "spaced_input.fa"] {
        let path = format!("../data/inputs/{file}");
        let data = read_sequence_file(Path::new(&path), b'%')?;
        for case in mask_cases(data.seq.clone(), data.start_positions.clone(), data.sequence_names.clone()) {
            assert_strategies_agree(&case, 1)?;
            assert_strategies_agree(&case, 3)?;
        }
    }
    Ok(())
}

// --------------------------------------------------
#[test]
fn radix_agrees_with_merge_on_random_dna() -> Result<()> {
    let mut rng = Rng(1);
    for long_n_run in [false, true] {
        let (text, starts, names) = random_genome(&mut rng, long_n_run);
        for case in mask_cases(text, starts, names) {
            assert_strategies_agree(&case, 4)?;
        }
    }
    Ok(())
}

// --------------------------------------------------
#[test]
fn radix_agrees_with_merge_on_protein() -> Result<()> {
    let mut rng = Rng(2);
    let alphabet = b"ACDEFGHIKLMNPQRSTVWY";
    let mut text = vec![];
    let mut starts = vec![];
    let mut names = vec![];
    for rec in 0..5 {
        if rec > 0 {
            text.push(b'%');
        }
        starts.push(text.len());
        names.push(format!("prot{rec}"));
        let len = between(&mut rng, 50, 400);
        text.extend((0..len).map(|_| alphabet[rng.below(alphabet.len())]));
    }
    text.push(b'$');

    // 22 symbols -> 5 bits each; weight 12 x 5 = 60 bits still fits
    for mask in ["101", "1101101", "111010011011", "1110110110111"] {
        let case = Case {
            text: text.clone(),
            starts: starts.clone(),
            names: names.clone(),
            is_dna: false,
            allow_ambiguity: false,
            seed_mask: Some(mask.to_string()),
            max_query_len: None,
        };
        assert_strategies_agree(&case, 4)?;
    }
    Ok(())
}

// --------------------------------------------------
/// In max-query-len mode the merge path's `find_lcp` can report LCPs
/// longer than the MQL after skipping known-equal characters, which
/// changes the order among suffixes with equal keys. So compare the radix
/// output with the independent definition, and check that merge agrees
/// with it up to tie order: same key at every rank.
#[test]
fn radix_matches_definition_on_max_query_len() -> Result<()> {
    let mut rng = Rng(3);
    let (text, starts, names) = random_genome(&mut rng, false);
    for mql in [1, 2, 5, 8, 12] {
        for (is_dna, allow_ambiguity) in [(true, false), (true, true), (false, false)] {
            let case = Case {
                text: text.clone(),
                starts: starts.clone(),
                names: names.clone(),
                is_dna,
                allow_ambiguity,
                seed_mask: None,
                max_query_len: Some(mql),
            };
            assert_matches_definition(&case, 4)?;

            let offsets = key_offsets(&case);
            let (merge_sa, _) = build(&case, SortStrategy::Merge, 1)?;
            let (radix_sa, _) = build(&case, SortStrategy::RadixPartitioned, 4)?;
            let keys = |sa: &[u32]| -> Vec<Vec<Option<u8>>> {
                sa.iter().map(|&p| masked(&text, &offsets, p as usize)).collect()
            };
            assert_eq!(keys(&radix_sa), keys(&merge_sa), "key order, {}", case.describe());
        }
    }
    Ok(())
}

// --------------------------------------------------
/// The merge path's suffix array should not depend on the number of
/// partitions. (In seed-mask mode `is_less` treats a suffix that matches a
/// pivot on all but the last care position as equal to the pivot, so such
/// suffixes all land in the pivot's partition; that skews partition sizes
/// but keeps the final order correct, which this test confirms.)
#[test]
fn merge_sa_is_independent_of_partition_count() -> Result<()> {
    let mut rng = Rng(9);
    let (text, starts, names) = random_genome(&mut rng, false);
    for mask in ["11011", "11101101101111"] {
        let case = Case {
            text: text.clone(),
            starts: starts.clone(),
            names: names.clone(),
            is_dna: true,
            allow_ambiguity: false,
            seed_mask: Some(mask.to_string()),
            max_query_len: None,
        };
        let (one, _) = build(&case, SortStrategy::Merge, 1)?;
        let (many, _) = build(&case, SortStrategy::Merge, 4)?;
        assert_eq!(many, one, "merge with 4 partitions vs 1, {}", case.describe());
    }
    Ok(())
}

// --------------------------------------------------
#[test]
fn radix_rejects_unsupported_modes() -> Result<()> {
    let mut rng = Rng(4);
    let (text, starts, names) = random_genome(&mut rng, true);

    // Full sort (no mask, no MQL)
    let case = Case {
        text: text.clone(),
        starts: starts.clone(),
        names: names.clone(),
        is_dna: true,
        allow_ambiguity: false,
        seed_mask: None,
        max_query_len: None,
    };
    assert!(build(&case, SortStrategy::RadixInMemory, 1).is_err());

    // MQL with long N runs and ambiguity allowed has special LCP semantics
    let case = Case {
        max_query_len: Some(5),
        allow_ambiguity: true,
        ..case.clone()
    };
    assert!(build(&case, SortStrategy::RadixPartitioned, 2).is_err());

    // Key too wide for 64 bits
    let case = Case {
        max_query_len: Some(30),
        allow_ambiguity: false,
        ..case.clone()
    };
    assert!(build(&case, SortStrategy::RadixInMemory, 1).is_err());
    Ok(())
}

// --------------------------------------------------
#[test]
fn radix_edge_cases() -> Result<()> {
    // Text shorter than the mask span; mask longer than the text
    let shorts: [&[u8]; 4] = [b"$", b"A$", b"AC$", b"ACGTN%AC$"];
    for text in shorts {
        for mask in ["101", "11000111", "1000000000001"] {
            let case = Case {
                text: text.to_vec(),
                starts: vec![0],
                names: vec!["s".to_string()],
                is_dna: true,
                allow_ambiguity: false,
                seed_mask: Some(mask.to_string()),
                max_query_len: None,
            };
            assert_matches_definition(&case, 2)?;
            assert_strategies_agree(&case, 2)?;
        }
    }

    // Repeat-heavy input: 300 near-copies of a 200 bp unit
    let mut rng = Rng(5);
    let unit = random_dna(&mut rng, 200);
    let mut text = vec![];
    for _ in 0..300 {
        let mut copy = unit.clone();
        if rng.below(4) == 0 {
            let at = rng.below(copy.len());
            copy[at] = b"ACGT"[rng.below(4)];
        }
        text.extend(copy);
    }
    text.push(b'$');
    for mask in ["11011", "11101101101111"] {
        let case = Case {
            text: text.clone(),
            starts: vec![0],
            names: vec!["rep".to_string()],
            is_dna: true,
            allow_ambiguity: false,
            seed_mask: Some(mask.to_string()),
            max_query_len: None,
        };
        assert_strategies_agree(&case, 4)?;
    }
    Ok(())
}

// --------------------------------------------------
/// Independent definition of the seed-mask order: compare the byte
/// sequences at the care offsets, with "past the end" sorting first, and
/// break ties by descending position.
#[test]
fn radix_output_is_sorted_and_complete() -> Result<()> {
    let mut rng = Rng(6);
    let (text, starts, names) = random_genome(&mut rng, true);
    for mask in ["101", "11011", "11101101101111"] {
        for (is_dna, allow_ambiguity) in [(true, false), (true, true), (false, false)] {
            let case = Case {
                text: text.clone(),
                starts: starts.clone(),
                names: names.clone(),
                is_dna,
                allow_ambiguity,
                seed_mask: Some(mask.to_string()),
                max_query_len: None,
            };
            assert_matches_definition(&case, 3)?;
        }
    }
    Ok(())
}

// --------------------------------------------------
#[test]
fn radix_index_searches_like_merge_index() -> Result<()> {
    let mut rng = Rng(8);
    let (text, starts, names) = random_genome(&mut rng, false);
    let mut queries = vec![];
    for _ in 0..40 {
        let len = between(&mut rng, 3, 16);
        let at = rng.below(text.len() - len - 1);
        queries.push(String::from_utf8(text[at..at + len].to_vec())?);
    }
    queries.push("ACGTACGTACGTACGT".to_string());

    for mask in ["11011", "11101101101111"] {
        let case = Case {
            text: text.clone(),
            starts: starts.clone(),
            names: names.clone(),
            is_dna: true,
            allow_ambiguity: false,
            seed_mask: Some(mask.to_string()),
            max_query_len: None,
        };
        let mut results = vec![];
        for strategy in [SortStrategy::Merge, SortStrategy::RadixPartitioned] {
            let outfile = NamedTempFile::new()?;
            let outpath = outfile.path().to_string_lossy().to_string();
            SufrBuilder::<u32>::new(case.args(outpath.clone(), 3, strategy))?;
            let mut sufr_file: SufrFile<u32> = SufrFile::read(&outpath, false)?;
            let counts = sufr_file.count(CountOptions {
                queries: queries.clone(),
                max_query_len: None,
                low_memory: true,
            })?;
            results.push(counts);
        }
        assert_eq!(results[0], results[1], "counts differ for mask {mask}");
        assert!(results[0].iter().any(|r| r.count > 0));
    }
    Ok(())
}
