use anyhow::Result;
use assert_cmd::Command;
use libsufr::{OpenOptions, Storage, SuffixArray, OUTFILE_VERSION};
use pretty_assertions::assert_eq;
use regex::Regex;
use std::fs;
use tempfile::{NamedTempFile, TempDir};

const PRG: &str = "sufr";
const SEQ1: &str = "../data/inputs/1.fa";
const SEQ2: &str = "../data/inputs/2.fa";
const SEQ3: &str = "../data/inputs/3.fa";
const LONG: &str = "../data/inputs/long_dna_sequence.fa";
const UNIPROT: &str = "../data/inputs/uniprot.fa";
const SUFR1: &str = "../data/expected/1.sufr";
const SUFR2: &str = "../data/expected/2.sufr";
const SUFR3: &str = "../data/expected/3.sufr";

struct ExtractOptions {
    queries: Vec<String>,
    prefix_len: Option<usize>,
    suffix_len: Option<usize>,
    max_query_len: Option<usize>,
}

struct ListOptions {
    show_rank: bool,
    show_suffix: bool,
    show_lcp: bool,
    suffix_len: Option<usize>,
    number: Option<usize>,
}

struct LocateOptions {
    queries: Vec<String>,
    absolute: bool,
    max_query_len: Option<usize>,
}

// --------------------------------------------------
#[test]
fn create_empty_input_fails() -> Result<()> {
    let outfile = NamedTempFile::new()?;
    let outpath = outfile.path().to_string_lossy().to_string();
    Command::cargo_bin(PRG)?
        .args(["create", "--dna", "-o", &outpath, "../data/inputs/empty.fa"])
        .assert()
        .failure();
    assert_eq!(
        fs::metadata(&outpath)?.len(),
        0,
        "nothing written for an empty input"
    );
    Ok(())
}

// --------------------------------------------------
fn count(
    filename: &str,
    queries: &[&str],
    expected_stdout: &str,
    expected_stderr: Option<&str>,
) -> Result<()> {
    for memory in &["", "-l", "-v"] {
        let mut args = vec!["count".to_string(), filename.to_string()];

        if !memory.is_empty() {
            args.push(memory.to_string())
        }

        for query in queries {
            args.push(query.to_string());
        }

        let output = Command::cargo_bin(PRG)?.args(&args).output().expect("fail");
        assert!(output.status.success());

        let stdout = String::from_utf8(output.stdout).expect("invalid UTF-8");
        assert_eq!(stdout, expected_stdout);

        if let Some(expected) = expected_stderr {
            let stderr = String::from_utf8(output.stderr).expect("invalid UTF-8");
            assert_eq!(stderr, expected);
        }
    }

    Ok(())
}

// --------------------------------------------------
#[test]
fn count_seq1() -> Result<()> {
    // cargo run -- co data/inputs/1.sufr AC X GT
    count(SUFR1, &["AC", "X", "GT"], "AC 2\nX 0\nGT 2\n", None)
}

// --------------------------------------------------
#[test]
fn count_seq3() -> Result<()> {
    // cargo run -- co data/inputs/3.sufr AAAAAAA TGTCTC TGATAGCAGCTTCTGAACTGGTTACCTGCCGT
    count(
        SUFR3,
        &["AAAAAAA", "TGTCTC", "TGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGT"],
        &[
            "AAAAAAA 1",
            "TGTCTC 1",
            "TGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGT 1",
            "",
        ]
        .join("\n"),
        None,
    )
}

// --------------------------------------------------
fn extract(
    filename: &str,
    opts: ExtractOptions,
    expected: &str,
    error: Option<&str>,
) -> Result<()> {
    // TODO: a search with -m writes a subsample cache under ~/.sufr;
    //       HOME points at a temp dir until sufr-22b9c3 decides
    //       where, or whether, that cache lives
    let home = TempDir::new()?;
    for memory in &["", "-l", "-v"] {
        let mut args = vec!["extract".to_string(), filename.to_string()];

        if !memory.is_empty() {
            args.push(memory.to_string());
        }

        if let Some(prefix_len) = opts.prefix_len {
            args.push("-p".to_string());
            args.push(prefix_len.to_string());
        }

        if let Some(suffix_len) = opts.suffix_len {
            args.push("-s".to_string());
            args.push(suffix_len.to_string());
        }

        if let Some(max_query_len) = opts.max_query_len {
            args.push("-m".to_string());
            args.push(max_query_len.to_string());
        }

        for query in &opts.queries {
            args.push(query.to_string());
        }

        let output = Command::cargo_bin(PRG)?
            .env("HOME", home.path())
            .args(&args)
            .output()
            .expect("fail");
        assert!(output.status.success());

        let stdout = String::from_utf8(output.stdout).expect("invalid UTF-8");
        assert_eq!(stdout, expected);

        if let Some(err) = error {
            let stderr = String::from_utf8(output.stderr).expect("invalid UTF-8");
            assert_eq!(stderr, err);
        }
    }

    Ok(())
}

// --------------------------------------------------
#[test]
fn extract_seq1() -> Result<()> {
    // cargo run -- ex data/inputs/1.sufr AC GT XX
    extract(
        SUFR1,
        ExtractOptions {
            queries: vec!["AC".to_string(), "GT".to_string(), "XX".to_string()],
            prefix_len: None,
            suffix_len: None,
            max_query_len: None,
        },
        &[
            ">1:6-10 AC 0",
            "ACGT",
            ">1:0-10 AC 0",
            "ACGTNNACGT",
            ">1:8-10 GT 0",
            "GT",
            ">1:2-10 GT 0",
            "GTNNACGT",
            "",
        ]
        .join("\n"),
        Some("XX not found\n"),
    )
}

// --------------------------------------------------
#[test]
fn extract_seq1_prefix_1() -> Result<()> {
    extract(
        SUFR1,
        ExtractOptions {
            queries: vec!["AC".to_string(), "GT".to_string()],
            prefix_len: Some(1),
            suffix_len: None,
            max_query_len: None,
        },
        &[
            ">1:5-10 AC 1",
            "NACGT",
            ">1:0-10 AC 0",
            "ACGTNNACGT",
            ">1:7-10 GT 1",
            "CGT",
            ">1:1-10 GT 1",
            "CGTNNACGT",
            "",
        ]
        .join("\n"),
        None,
    )
}

// --------------------------------------------------
#[test]
fn extract_seq1_suf_3() -> Result<()> {
    extract(
        SUFR1,
        ExtractOptions {
            queries: vec!["AC".to_string(), "GT".to_string()],
            prefix_len: None,
            suffix_len: Some(3),
            max_query_len: None,
        },
        &[
            ">1:6-9 AC 0",
            "ACG",
            ">1:0-3 AC 0",
            "ACG",
            ">1:8-10 GT 0",
            "GT",
            ">1:2-5 GT 0",
            "GTN",
            "",
        ]
        .join("\n"),
        None,
    )
}

// --------------------------------------------------
#[test]
fn extract_seq1_pre_1_suf_3() -> Result<()> {
    extract(
        SUFR1,
        ExtractOptions {
            queries: vec!["AC".to_string(), "GT".to_string()],
            prefix_len: Some(1),
            suffix_len: Some(3),
            max_query_len: None,
        },
        &[
            ">1:5-9 AC 1",
            "NACG",
            ">1:0-3 AC 0",
            "ACG",
            ">1:7-10 GT 1",
            "CGT",
            ">1:1-5 GT 1",
            "CGTN",
            "",
        ]
        .join("\n"),
        None,
    )
}

// --------------------------------------------------
#[test]
fn extract_uniprot() -> Result<()> {
    extract(
        "../data/expected/uniprot.sufr",
        ExtractOptions {
            queries: vec![
                "RNELNNEEA".to_string(),
                "DTPTNCPT".to_string(),
                "GSGLSLLSD".to_string(),
            ],
            prefix_len: Some(4),
            suffix_len: Some(12),
            max_query_len: None,
        },
        &[
            ">sp|Q9U408|14331_ECHGR:38-54 RNELNNEEA 4",
            "MAKMRNELNNEEANLL",
            ">sp|Q6GZX3|002L_FRG3G:218-234 DTPTNCPT 4",
            "GTQRDTPTNCPTQVCQ",
            ">sp|Q6GZW6|009L_FRG3G:390-406 GSGLSLLSD 4",
            "EYVNGSGLSLLSDILL",
            "",
        ]
        .join("\n"),
        None,
    )
}

// --------------------------------------------------
#[test]
fn extract_uniprot_mql() -> Result<()> {
    // In this test, the MQL 3 means all the found suffixes start with "RNE"
    // cargo run -- ex data/expected/uniprot.sufr RNELNNEEA -s 10 -m 3
    extract(
        "../data/expected/uniprot.sufr",
        ExtractOptions {
            queries: vec!["RNELNNEEA".to_string()],
            prefix_len: None,
            suffix_len: Some(10),
            max_query_len: Some(3),
        },
        &[
            ">sp|Q6GZW1|014R_FRG3G:111-119 RNELNNEEA 0",
            "RNEEDDDG",
            ">sp|P0C9G4|1101L_ASFP4:80-90 RNELNNEEA 0",
            "RNEFCTYYVT",
            ">sp|P0C9G1|1101L_ASFWA:80-90 RNELNNEEA 0",
            "RNEFCTYYVT",
            ">sp|Q9U408|14331_ECHGR:42-52 RNELNNEEA 0",
            "RNELNNEEAN",
            ">sp|O55726|110R_IIV6:31-41 RNELNNEEA 0",
            "RNEPSHYQTV",
            ">sp|Q196U3|117L_IIV3:27-37 RNELNNEEA 0",
            "RNEYDNAVAS",
            ">sp|Q197E9|011L_IIV3:55-65 RNELNNEEA 0",
            "RNEYNKVHIE",
            "",
        ]
        .join("\n"),
        None,
    )
}

// --------------------------------------------------
#[test]
fn extract_uniprot_masked() -> Result<()> {
    // In this test, the MQL 3 means all the found suffixes start with "R*EL"
    // cargo run -- ex -s 10 -m 3 data/expected/uniprot-masked.sufr RNELNNEEA
    extract(
        "../data/expected/uniprot-masked.sufr",
        ExtractOptions {
            queries: vec!["RNELNNEEA".to_string()],
            prefix_len: None,
            suffix_len: Some(10),
            max_query_len: Some(3),
        },
        &[
            ">sp|P32234|128UP_DROME:38-48 RNELNNEEA 0",
            "RRELISPKGG",
            ">sp|P26709|1107L_ASFL5:128-137 RNELNNEEA 0",
            "RKELKKDEF",
            ">sp|P0DO85|10H_STRNX:151-161 RNELNNEEA 0",
            "RLELLKHIRV",
            ">sp|Q9U408|14331_ECHGR:42-52 RNELNNEEA 0",
            "RNELNNEEAN",
            ">sp|P19084|11S3_HELAN:355-365 RNELNNEEA 0",
            "RGELRPNAIQ",
            ">sp|Q6GZV8|017L_FRG3G:18-28 RNELNNEEA 0",
            "RGELSALSAA",
            ">sp|Q6GZW6|009L_FRG3G:653-663 RNELNNEEA 0",
            "RLELSAPYGS",
            "",
        ]
        .join("\n"),
        None,
    )
}

// --------------------------------------------------
fn list(filename: &str, opts: ListOptions, expected: &str) -> Result<()> {
    for memory in &["", "-v"] {
        let mut args = vec!["list".to_string(), filename.to_string()];

        if !memory.is_empty() {
            args.push(memory.to_string());
        }

        if opts.show_suffix {
            args.push("-s".to_string());
        }

        if opts.show_rank {
            args.push("-r".to_string());
        }

        if opts.show_lcp {
            args.push("-p".to_string());
        }

        if let Some(len) = opts.suffix_len {
            args.push("--len".to_string());
            args.push(len.to_string());
        }

        if let Some(num) = opts.number {
            args.push("-n".to_string());
            args.push(num.to_string());
        }

        let output = Command::cargo_bin(PRG)?.args(&args).output().expect("fail");
        assert!(output.status.success());

        let stdout = String::from_utf8(output.stdout).expect("invalid UTF-8");
        assert_eq!(stdout, expected);
    }

    Ok(())
}

// --------------------------------------------------
#[test]
fn list_sufr1_no_opts() -> Result<()> {
    list(
        SUFR1,
        ListOptions {
            show_lcp: false,
            show_suffix: false,
            show_rank: false,
            suffix_len: None,
            number: None,
        },
        &[
            "$",
            "ACGT$",
            "ACGTNNACGT$",
            "CGT$",
            "CGTNNACGT$",
            "GT$",
            "GTNNACGT$",
            "T$",
            "TNNACGT$",
            "",
        ]
        .join("\n"),
    )
}

// --------------------------------------------------
#[test]
fn list_sufr1_show_suffix() -> Result<()> {
    list(
        SUFR1,
        ListOptions {
            show_lcp: false,
            show_suffix: true,
            show_rank: false,
            suffix_len: None,
            number: None,
        },
        &[
            "10 $",
            " 6 ACGT$",
            " 0 ACGTNNACGT$",
            " 7 CGT$",
            " 1 CGTNNACGT$",
            " 8 GT$",
            " 2 GTNNACGT$",
            " 9 T$",
            " 3 TNNACGT$",
            "",
        ]
        .join("\n"),
    )
}

// --------------------------------------------------
#[test]
fn list_sufr1_show_lcp() -> Result<()> {
    list(
        SUFR1,
        ListOptions {
            show_lcp: true,
            show_suffix: false,
            show_rank: false,
            suffix_len: None,
            number: None,
        },
        &[
            " 0 $",
            " 0 ACGT$",
            " 4 ACGTNNACGT$",
            " 0 CGT$",
            " 3 CGTNNACGT$",
            " 0 GT$",
            " 2 GTNNACGT$",
            " 0 T$",
            " 1 TNNACGT$",
            "",
        ]
        .join("\n"),
    )
}

// --------------------------------------------------
#[test]
fn list_sufr1_show_rank_suffix_lcp() -> Result<()> {
    // cargo run -- ls -Lsr data/expected/1.sufr
    list(
        SUFR1,
        ListOptions {
            show_lcp: true,
            show_suffix: true,
            show_rank: true,
            suffix_len: None,
            number: None,
        },
        &[
            " 0 10  0 $",
            " 1  6  0 ACGT$",
            " 2  0  4 ACGTNNACGT$",
            " 3  7  0 CGT$",
            " 4  1  3 CGTNNACGT$",
            " 5  8  0 GT$",
            " 6  2  2 GTNNACGT$",
            " 7  9  0 T$",
            " 8  3  1 TNNACGT$",
            "",
        ]
        .join("\n"),
    )
}

// --------------------------------------------------
#[test]
fn list_sufr1_show_rank() -> Result<()> {
    list(
        SUFR1,
        ListOptions {
            show_lcp: false,
            show_suffix: false,
            show_rank: true,
            suffix_len: None,
            number: None,
        },
        &[
            " 0 $",
            " 1 ACGT$",
            " 2 ACGTNNACGT$",
            " 3 CGT$",
            " 4 CGTNNACGT$",
            " 5 GT$",
            " 6 GTNNACGT$",
            " 7 T$",
            " 8 TNNACGT$",
            "",
        ]
        .join("\n"),
    )
}

// --------------------------------------------------
#[test]
fn list_sufr1_len_3() -> Result<()> {
    // cargo run -- ls -l 3 data/expected/1.sufr
    list(
        SUFR1,
        ListOptions {
            show_lcp: false,
            show_suffix: false,
            show_rank: false,
            suffix_len: Some(3),
            number: None,
        },
        &[
            "$", "ACG", "ACG", "CGT", "CGT", "GT$", "GTN", "T$", "TNN", "",
        ]
        .join("\n"),
    )
}

// --------------------------------------------------
#[test]
fn list_sufr1_limit_3() -> Result<()> {
    // cargo run -- ls -n 3 data/expected/1.sufr
    list(
        SUFR1,
        ListOptions {
            show_lcp: false,
            show_suffix: false,
            show_rank: false,
            suffix_len: None,
            number: Some(3),
        },
        &["$", "ACGT$", "ACGTNNACGT$", ""].join("\n"),
    )
}

// --------------------------------------------------
fn locate(filename: &str, opts: LocateOptions, expected_file: &str) -> Result<()> {
    // TODO: as in extract, HOME points at a temp dir until sufr-22b9c3
    let home = TempDir::new()?;
    for memory in ["", "-l", "-v"] {
        let mut args = vec!["locate".to_string(), filename.to_string()];

        if opts.absolute {
            args.push("-a".to_string());
        }

        if !memory.is_empty() {
            args.push(memory.to_string());
        }

        if let Some(max_query_len) = opts.max_query_len {
            args.push("-m".to_string());
            args.push(max_query_len.to_string());
        }

        for query in &opts.queries {
            args.push(query.to_string());
        }

        let output = Command::cargo_bin(PRG)?
            .env("HOME", home.path())
            .args(&args)
            .output()
            .expect("fail");
        assert!(output.status.success());

        let actual = String::from_utf8(output.stdout).expect("invalid UTF-8");
        let expected = fs::read_to_string(expected_file)?;
        assert_eq!(actual, expected);
    }

    Ok(())
}

// --------------------------------------------------
#[test]
fn locate_seq1_relative() -> Result<()> {
    // cargo run -- lo data/expected/2.sufr AC GT
    // AC
    // ABC 0,4
    // DEF 0,4
    // //
    // GT
    // ABC 2,6
    // DEF 2,6
    // //

    locate(
        SUFR2,
        LocateOptions {
            queries: vec!["AC".to_string(), "GT".to_string()],
            absolute: false,
            max_query_len: None,
        },
        "../data/expected/locate1.out",
    )
}

// --------------------------------------------------
#[test]
fn locate_seq1_absolute() -> Result<()> {
    // cargo run -- lo data/expected/2.sufr AC GT -a
    // AC 13 4 9 0
    // GT 15 6 11 2

    locate(
        SUFR2,
        LocateOptions {
            queries: vec!["AC".to_string(), "GT".to_string()],
            absolute: true,
            max_query_len: None,
        },
        "../data/expected/locate-abs.out",
    )
}

// --------------------------------------------------
#[test]
fn locate_uniprot() -> Result<()> {
    // cargo run -- lo data/expected/uniprot.sufr RNELNNEEA DTPTNCPT GSGLSLLSD
    // RNELNNEEA
    // sp|Q9U408|14331_ECHGR 42
    // //
    // DTPTNCPT
    // sp|Q6GZX3|002L_FRG3G 222
    // //
    // GSGLSLLSD
    // sp|Q6GZW6|009L_FRG3G 394
    // //

    locate(
        "../data/expected/uniprot.sufr",
        LocateOptions {
            queries: vec![
                "RNELNNEEA".to_string(),
                "DTPTNCPT".to_string(),
                "GSGLSLLSD".to_string(),
            ],
            absolute: false,
            max_query_len: None,
        },
        "../data/expected/uniprot-search1.out",
    )
}

// --------------------------------------------------
#[test]
fn locate_uniprot_masked() -> Result<()> {
    // cargo run -- lo data/expected/uniprot-masked.sufr RNEL DTPT GSGL
    locate(
        "../data/expected/uniprot-masked.sufr",
        LocateOptions {
            queries: vec!["RNEL".to_string(), "DTPT".to_string(), "GSGL".to_string()],
            absolute: false,
            max_query_len: None,
        },
        "../data/expected/uniprot-search-masked.out",
    )
}

// --------------------------------------------------
#[test]
fn locate_uniprot_masked_max_query_len() -> Result<()> {
    // cargo run -- lo data/expected/uniprot-masked.sufr -m 3 RNELNNEEA DTPTNCPT GSGLSLLSD

    locate(
        "../data/expected/uniprot-masked.sufr",
        LocateOptions {
            queries: vec![
                "RNELNNEEA".to_string(),
                "DTPTNCPT".to_string(),
                "GSGLSLLSD".to_string(),
            ],
            absolute: false,
            max_query_len: Some(3),
        },
        "../data/expected/uniprot-search-masked-mql-3.out",
    )
}

// --------------------------------------------------
#[test]
fn locate_uniprot_masked_absolute() -> Result<()> {
    // cargo run -- lo data/expected/uniprot-masked.sufr RNEL -a
    // RNEL 54791 46515 37970 62005 52278 6386 4124

    locate(
        "../data/expected/uniprot-masked.sufr",
        LocateOptions {
            queries: vec!["RNEL".to_string()],
            absolute: true,
            max_query_len: None,
        },
        "../data/expected/uniprot-search-masked-absolute.out",
    )
}

// --------------------------------------------------
#[test]
fn locate_long_dna_low_memory() -> Result<()> {
    // cargo run -- lo data/expected/long_dna_sequence.sufr \
    // CATGTTGTCACG CCATGGGAC GGATGAAGAAAAGCA
    // CATGTTGTCACG
    // Seq1 2566
    // //
    // CCATGGGAC
    // Seq1 3014
    // //
    // GGATGAAGAAAAGCA
    // Seq1 1026
    // //

    locate(
        "../data/expected/long_dna_sequence.sufr",
        LocateOptions {
            queries: vec![
                "CATGTTGTCACG".to_string(),
                "CCATGGGAC".to_string(),
                "GGATGAAGAAAAGCA".to_string(),
            ],
            absolute: false,
            max_query_len: None,
        },
        "../data/expected/locate_long_dna.out",
    )
}

// --------------------------------------------------
#[test]
fn locate_long_dna_max_query_len() -> Result<()> {
    // cargo run -- lo data/expected/long_dna_sequence.sufr \
    // CATGTTGTCACG CCATGGGAC GGATGAAGAAAAGCA -m 6
    //
    // CATGTTGTCACG
    // Seq1 2566,13056,20444
    // //
    // CCATGGGAC
    // Seq1 3014
    // //
    // GGATGAAGAAAAGCA
    // Seq1 1026,2905,13253,13508,14250,14624,20465
    // //

    locate(
        "../data/expected/long_dna_sequence.sufr",
        LocateOptions {
            queries: vec![
                "CATGTTGTCACG".to_string(),
                "CCATGGGAC".to_string(),
                "GGATGAAGAAAAGCA".to_string(),
            ],
            absolute: false,
            max_query_len: Some(6),
        },
        "../data/expected/locate_long_dna_mql_6.out",
    )
}

// --------------------------------------------------
fn summarize(filename: &str, expected: Vec<(&str, &str)>) -> Result<()> {
    let args = vec!["summarize".to_string(), filename.to_string()];
    let output = Command::cargo_bin(PRG)?.args(&args).output().expect("fail");
    assert!(output.status.success());

    let stdout = String::from_utf8(output.stdout).expect("invalid UTF-8");

    for (row_name, value) in expected {
        let pattern = format!(r#"[|] {}\s+[|] ([^|]+)"#, row_name);
        let regex = Regex::new(&pattern).unwrap();
        let caps = regex
            .captures(&stdout)
            .unwrap_or_else(|| panic!("no row {row_name} in:\n{stdout}"));
        assert_eq!(
            caps.get(1).unwrap().as_str().trim(),
            value,
            "row {row_name}"
        );
    }

    Ok(())
}

// --------------------------------------------------
#[test]
fn summarize_sufr1() -> Result<()> {
    summarize(
        SUFR1,
        vec![
            ("File Size", "172 bytes"),
            ("File Version", &OUTFILE_VERSION.to_string()),
            ("DNA", "true"),
            ("Allow Ambiguity", "false"),
            ("Ignore Softmask", "false"),
            ("Text Length", "11"),
            ("Len Suffixes", "9"),
            ("Max query len", "0"),
            ("Num sequences", "1"),
            ("Sequence starts", "0"),
            ("Sequence names", "1"),
        ],
    )
}

// --------------------------------------------------
fn run_stdout(args: &[&str]) -> Result<(String, String)> {
    let output = Command::cargo_bin(PRG)?.args(args).output().expect("fail");
    assert!(output.status.success());
    Ok((
        String::from_utf8(output.stdout).expect("invalid UTF-8"),
        String::from_utf8(output.stderr).expect("invalid UTF-8"),
    ))
}

// --------------------------------------------------
#[test]
fn create_bytes_match_fixtures() -> Result<()> {
    let cases: Vec<(&str, &[&str], &str)> = vec![
        ("1", &["--dna"], SEQ1),
        ("2", &["--dna"], SEQ2),
        ("3", &["--dna"], SEQ3),
        ("2d", &["--dna", "--sequence-delimiter", "N"], SEQ2),
        ("abba", &[], "../data/inputs/abba.fa"),
        ("1n", &["--dna", "--allow-ambiguity"], SEQ1),
        ("2n", &["--dna", "--allow-ambiguity"], SEQ2),
        ("2s", &["--dna", "--ignore-softmask"], SEQ2),
        (
            "2ns",
            &["--dna", "--allow-ambiguity", "--ignore-softmask"],
            SEQ2,
        ),
        ("long_dna_sequence", &["--dna"], LONG),
        (
            "long_dna_sequence_allow_ambiguity",
            &["--dna", "--allow-ambiguity"],
            LONG,
        ),
        ("uniprot", &[], UNIPROT),
        ("uniprot-masked", &["-s", "10111011"], UNIPROT),
    ];

    for (name, flags, input) in cases {
        let outfile = NamedTempFile::new()?;
        let outpath = outfile.path().to_string_lossy().to_string();
        let mut args = vec!["create", "--write-lcp", "-o", &outpath, input];
        args.extend(flags);
        run_stdout(&args)?;

        let actual = fs::read(&outpath)?;
        let expected = fs::read(format!("../data/expected/{name}.sufr"))?;
        assert!(actual == expected, "bytes differ from fixture {name}");
    }

    Ok(())
}

// --------------------------------------------------
#[test]
fn create_num_partitions_same_bytes() -> Result<()> {
    let mut outputs = vec![];
    for parts in ["3", "1024"] {
        let outfile = NamedTempFile::new()?;
        let outpath = outfile.path().to_string_lossy().to_string();
        run_stdout(&[
            "create",
            "--dna",
            "--write-lcp",
            "-n",
            parts,
            "-o",
            &outpath,
            SEQ1,
        ])?;
        outputs.push(fs::read(&outpath)?);
    }
    assert!(outputs[0] == outputs[1]);
    assert!(outputs[0] == fs::read(SUFR1)?);
    Ok(())
}

// --------------------------------------------------
#[test]
fn create_max_query_len() -> Result<()> {
    let outfile = NamedTempFile::new()?;
    let outpath = outfile.path().to_string_lossy().to_string();
    run_stdout(&["create", "--dna", "-m", "3", "-o", &outpath, SEQ1])?;

    let actual = SuffixArray::open_with(
        &outpath,
        &OpenOptions {
            text: Storage::Memory,
            suffixes: Storage::Disk,
            max_query_len: None,
        },
    )?;
    let actual_sa = actual.suffixes(0..actual.header().num_suffixes)?;
    assert_eq!(actual_sa, [10, 6, 0, 7, 1, 8, 2, 9, 3]);
    Ok(())
}

// --------------------------------------------------
#[test]
fn create_default_output_name() -> Result<()> {
    let dir = TempDir::new()?;
    let input = fs::canonicalize(SEQ1)?;
    let output = Command::cargo_bin(PRG)?
        .current_dir(dir.path())
        .args(["create", "--dna", "--write-lcp"])
        .arg(&input)
        .output()
        .expect("fail");
    assert!(output.status.success());
    assert!(output.stdout.is_empty());

    let names: Vec<_> = fs::read_dir(dir.path())?
        .map(|e| e.map(|e| e.file_name().to_string_lossy().to_string()))
        .collect::<std::io::Result<_>>()?;
    assert_eq!(names, ["1.sufr"]);
    assert!(fs::read(dir.path().join("1.sufr"))? == fs::read(SUFR1)?);
    Ok(())
}

// --------------------------------------------------
fn check_output_file(subcommand: &str, filename: &str, queries: &[&str]) -> Result<()> {
    let mut args = vec![subcommand, filename];
    args.extend(queries);
    let (expected, _) = run_stdout(&args)?;
    assert!(!expected.is_empty());

    let outfile = NamedTempFile::new()?;
    let outpath = outfile.path().to_string_lossy().to_string();
    let mut args = vec![subcommand, "-o", &outpath, filename];
    args.extend(queries);
    let (stdout, _) = run_stdout(&args)?;

    assert_eq!(stdout, "");
    assert_eq!(fs::read_to_string(&outpath)?, expected);
    Ok(())
}

// --------------------------------------------------
#[test]
fn count_output_file() -> Result<()> {
    check_output_file("count", SUFR1, &["AC", "X", "GT"])
}

// --------------------------------------------------
#[test]
fn extract_output_file() -> Result<()> {
    check_output_file("extract", SUFR1, &["AC", "GT"])
}

// --------------------------------------------------
#[test]
fn locate_output_file() -> Result<()> {
    check_output_file("locate", SUFR2, &["AC", "GT"])
}

// --------------------------------------------------
#[test]
fn locate_not_found_goes_to_stderr() -> Result<()> {
    let outfile = NamedTempFile::new()?;
    let outpath = outfile.path().to_string_lossy().to_string();
    let (stdout, stderr) = run_stdout(&["locate", "-o", &outpath, SUFR1, "AC", "X"])?;
    assert_eq!(stdout, "");
    assert_eq!(stderr, "X not found\n");
    assert_eq!(fs::read_to_string(&outpath)?, "AC\n1 0,6\n//\n");
    Ok(())
}

// --------------------------------------------------
#[test]
fn list_ranks() -> Result<()> {
    let (stdout, _) = run_stdout(&["list", SUFR1, "0", "2-3"])?;
    assert_eq!(stdout, "$\nACGTNNACGT$\nCGT$\n");
    Ok(())
}

// --------------------------------------------------
#[test]
fn list_ranks_comma_with_columns() -> Result<()> {
    let (stdout, stderr) = run_stdout(&["list", "-r", "-s", "-p", SUFR1, "5-6,9"])?;
    assert_eq!(stdout, " 5  8  0 GT$\n 6  2  2 GTNNACGT$\n");
    assert_eq!(stderr, "Invalid rank: 9\n");
    Ok(())
}

// --------------------------------------------------
#[test]
fn list_ranks_len() -> Result<()> {
    let (stdout, _) = run_stdout(&["list", "--len", "3", SUFR1, "0", "2-3"])?;
    assert_eq!(stdout, "$\nACG\nCGT\n");
    Ok(())
}

// --------------------------------------------------
#[test]
fn list_short_l_is_rejected() -> Result<()> {
    Command::cargo_bin(PRG)?
        .args(["list", "-l", SUFR1, "0"])
        .assert()
        .failure();
    Ok(())
}

// --------------------------------------------------
#[test]
fn list_out_of_range_rank_warns_and_succeeds() -> Result<()> {
    let (stdout, stderr) = run_stdout(&["list", SUFR1, "99"])?;
    assert_eq!(stdout, "");
    assert_eq!(stderr, "Invalid rank: 99\n");
    Ok(())
}

// --------------------------------------------------
#[test]
fn list_output_flag_writes_the_file() -> Result<()> {
    let (expected, _) = run_stdout(&["list", SUFR1, "0", "1"])?;

    let outfile = NamedTempFile::new()?;
    let outpath = outfile.path().to_string_lossy().to_string();
    let (stdout, _) = run_stdout(&["list", "-o", &outpath, SUFR1, "0", "1"])?;

    assert_eq!(stdout, "");
    assert_eq!(fs::read_to_string(&outpath)?, expected);
    Ok(())
}

// --------------------------------------------------
#[test]
fn count_query_file() -> Result<()> {
    let queries = NamedTempFile::new()?;
    fs::write(queries.path(), "AC GT\nX\n")?;
    let path = queries.path().to_string_lossy().to_string();
    let (stdout, _) = run_stdout(&["count", SUFR1, &path])?;
    assert_eq!(stdout, "AC 2\nGT 2\nX 0\n");

    let (stdout, _) = run_stdout(&["count", SUFR1, "GT", &path])?;
    assert_eq!(stdout, "GT 2\nAC 2\nGT 2\nX 0\n");
    Ok(())
}

// --------------------------------------------------
// a cap changes which suffixes a search can tell apart, never how
// many match a query shorter than it
#[test]
fn count_under_a_cap_matches_count_without_one() -> Result<()> {
    // TODO: as in extract, HOME points at a temp dir until sufr-22b9c3
    let home = TempDir::new()?;
    let outfile = NamedTempFile::new()?;
    let outpath = outfile.path().to_string_lossy().to_string();
    Command::cargo_bin(PRG)?
        .args(["create", "--dna", "-o", &outpath, "../data/inputs/smol.fa"])
        .assert()
        .success();

    for cap in ["", "-m 2", "-m 3"] {
        let mut args = vec!["count", &outpath, "A", "AC"];
        args.extend(cap.split_whitespace());
        let output = Command::cargo_bin(PRG)?
            .env("HOME", home.path())
            .args(&args)
            .output()?;
        assert!(output.status.success());
        assert_eq!(String::from_utf8(output.stdout)?, "A 93\nAC 11\n", "{cap}");
    }

    // the T group is the last one under -m 1
    let output = Command::cargo_bin(PRG)?
        .env("HOME", home.path())
        .args(["count", "-m", "1", SUFR1, "T"])
        .output()?;
    assert!(output.status.success());
    assert_eq!(String::from_utf8(output.stdout)?, "T 2\n");
    Ok(())
}
// --------------------------------------------------
#[test]
fn locate_query_file() -> Result<()> {
    let queries = NamedTempFile::new()?;
    fs::write(queries.path(), "AC GT\nX\n")?;
    let path = queries.path().to_string_lossy().to_string();
    let (stdout, stderr) = run_stdout(&["locate", SUFR1, &path])?;
    assert_eq!(stdout, "AC\n1 0,6\n//\nGT\n1 2,8\n//\n");
    assert_eq!(stderr, "X not found\n");
    Ok(())
}
