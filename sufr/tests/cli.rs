use anyhow::Result;
use assert_cmd::Command;
use libsufr::{
    sufr_file::SufrFile,
    types::{SeedMask, OUTFILE_VERSION},
};
use pretty_assertions::assert_eq;
use regex::Regex;
use std::{
    fs::{self, File},
    io::{BufRead, BufReader},
};
use tempfile::NamedTempFile;

const PRG: &str = "sufr";
const SEQ1: &str = "../data/inputs/1.fa";
const SEQ2: &str = "../data/inputs/2.fa";
const SEQ3: &str = "../data/inputs/3.fa";
const LONG: &str = "../data/inputs/long_dna_sequence.fa";
const UNIPROT: &str = "../data/inputs/uniprot.fa";
const SUFR1: &str = "../data/expected/1.sufr";
const SUFR2: &str = "../data/expected/2.sufr";
const SUFR3: &str = "../data/expected/3.sufr";

struct CreateOptions {
    is_dna: bool,
    allow_ambiguity: bool,
    ignore_softmask: bool,
    sequence_delimiter: Option<char>,
    seed_mask: Option<String>,
}

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
    low_memory: bool,
    max_query_len: Option<usize>,
}

// --------------------------------------------------
fn create(input_file: &str, expected_file: &str, opts: CreateOptions) -> Result<()> {
    let outfile = NamedTempFile::new()?;
    let outpath = &outfile.path().to_string_lossy();
    let mut args = vec![
        "create".to_string(),
        "-o".to_string(),
        outpath.to_string(),
        input_file.to_string(),
    ];

    if opts.is_dna {
        args.push("--dna".to_string());
    }

    if opts.allow_ambiguity {
        args.push("--allow-ambiguity".to_string());
    }

    if opts.ignore_softmask {
        args.push("--ignore-softmask".to_string());
    }

    if let Some(delim) = opts.sequence_delimiter {
        let mut tmp = vec!["--sequence-delimiter".to_string(), delim.to_string()];
        args.append(&mut tmp);
    }

    if let Some(mask) = opts.seed_mask {
        let mut tmp = vec!["--seed-mask".to_string(), mask.to_string()];
        args.append(&mut tmp);
    }

    let output = Command::cargo_bin(PRG)?.args(args).output().expect("fail");

    assert!(output.status.success());
    assert!(outfile.path().exists());

    let mut actual: SufrFile<u32> = SufrFile::read(outpath, false)?;
    let mut expected: SufrFile<u32> = SufrFile::read(expected_file, false)?;

    let actual_sa: Vec<_> = actual.suffix_array_file.iter().collect();
    let expected_sa: Vec<_> = expected.suffix_array_file.iter().collect();
    assert_eq!(actual_sa, expected_sa);

    Ok(())
}

// --------------------------------------------------
// TODO: revive once check is useful
//fn check(filename: &str) -> Result<()> {
//    let output = Command::cargo_bin(PRG)?
//        .args(["check", filename])
//        .output()
//        .expect("fail");
//
//    assert!(output.status.success());
//
//    let stdout = String::from_utf8(output.stdout).expect("invalid UTF-8");
//
//    assert!(stdout.contains("found 0 errors"));
//
//    Ok(())
//}
//
// --------------------------------------------------
//#[test]
//fn check_seq1() -> Result<()> {
//    check("../data/expected/1.sufr")
//}

// --------------------------------------------------
#[test]
fn create_empty_dies() -> Result<()> {
    Command::cargo_bin(PRG)?
        .args(["create", "../data/expected/empty.sa"])
        .assert()
        .failure();
    Ok(())
}

// --------------------------------------------------
#[test]
fn create_seq1() -> Result<()> {
    create(
        SEQ1,
        "../data/expected/1.sufr",
        CreateOptions {
            is_dna: true,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_delimiter: None,
            seed_mask: None,
        },
    )
}

// --------------------------------------------------
#[test]
fn create_seq1_allow_ambiguity() -> Result<()> {
    create(
        SEQ1,
        "../data/expected/1n.sufr",
        CreateOptions {
            is_dna: true,
            allow_ambiguity: true,
            ignore_softmask: false,
            sequence_delimiter: None,
            seed_mask: None,
        },
    )
}

// --------------------------------------------------
#[test]
fn create_seq2() -> Result<()> {
    create(
        SEQ2,
        "../data/expected/2.sufr",
        CreateOptions {
            is_dna: true,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_delimiter: None,
            seed_mask: None,
        },
    )
}

// --------------------------------------------------
#[test]
fn create_seq2_sequence_delimiter() -> Result<()> {
    create(
        SEQ2,
        "../data/expected/2d.sufr",
        CreateOptions {
            is_dna: true,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_delimiter: Some('N'),
            seed_mask: None,
        },
    )
}

// --------------------------------------------------
#[test]
fn create_seq2_allow_ambiguity() -> Result<()> {
    create(
        SEQ2,
        "../data/expected/2n.sufr",
        CreateOptions {
            is_dna: true,
            allow_ambiguity: true,
            ignore_softmask: false,
            sequence_delimiter: None,
            seed_mask: None,
        },
    )
}

// --------------------------------------------------
#[test]
fn create_seq2_ignore_softmask() -> Result<()> {
    create(
        SEQ2,
        "../data/expected/2s.sufr",
        CreateOptions {
            is_dna: true,
            allow_ambiguity: false,
            ignore_softmask: true,
            sequence_delimiter: None,
            seed_mask: None,
        },
    )
}

// --------------------------------------------------
#[test]
fn create_seq2_allow_ambiguity_ignore_softmask() -> Result<()> {
    create(
        SEQ2,
        "../data/expected/2ns.sufr",
        CreateOptions {
            is_dna: true,
            allow_ambiguity: true,
            ignore_softmask: true,
            sequence_delimiter: None,
            seed_mask: None,
        },
    )
}

// --------------------------------------------------
#[test]
fn create_seq3() -> Result<()> {
    create(
        SEQ3,
        "../data/expected/3.sufr",
        CreateOptions {
            is_dna: true,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_delimiter: None,
            seed_mask: None,
        },
    )
}

// --------------------------------------------------
#[test]
fn create_long_dna() -> Result<()> {
    create(
        LONG,
        "../data/expected/long_dna_sequence.sufr",
        CreateOptions {
            is_dna: true,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_delimiter: None,
            seed_mask: None,
        },
    )
}

// --------------------------------------------------
#[test]
fn create_long_dna_allow_ambiguity() -> Result<()> {
    create(
        LONG,
        "../data/expected/long_dna_sequence_allow_ambiguity.sufr",
        CreateOptions {
            is_dna: true,
            allow_ambiguity: true,
            ignore_softmask: false,
            sequence_delimiter: None,
            seed_mask: None,
        },
    )
}

// --------------------------------------------------
#[test]
fn create_protein() -> Result<()> {
    create(
        UNIPROT,
        "../data/expected/uniprot.sufr",
        CreateOptions {
            is_dna: false,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_delimiter: None,
            seed_mask: None,
        },
    )
}

// --------------------------------------------------
#[test]
fn create_protein_masked() -> Result<()> {
    create(
        UNIPROT,
        "../data/expected/uniprot-masked.sufr",
        CreateOptions {
            is_dna: false,
            allow_ambiguity: false,
            ignore_softmask: false,
            sequence_delimiter: None,
            seed_mask: Some("10111011".to_string()),
        },
    )
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

        dbg!(&args);
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

        let output = Command::cargo_bin(PRG)?.args(&args).output().expect("fail");
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
            ">1:6-11 AC 0",
            "ACGT$",
            ">1:0-11 AC 0",
            "ACGTNNACGT$",
            ">1:8-11 GT 0",
            "GT$",
            ">1:2-11 GT 0",
            "GTNNACGT$",
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
            ">1:5-11 AC 1",
            "NACGT$",
            ">1:0-11 AC 0",
            "ACGTNNACGT$",
            ">1:7-11 GT 1",
            "CGT$",
            ">1:1-11 GT 1",
            "CGTNNACGT$",
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
            ">1:8-11 GT 0",
            "GT$",
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
            ">1:7-11 GT 1",
            "CGT$",
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
            ">sp|Q6GZW1|014R_FRG3G:111-121 RNELNNEEA 0",
            "RNEEDDDG%M",
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
            ">sp|P26709|1107L_ASFL5:128-138 RNELNNEEA 0",
            "RKELKKDEF%",
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
#[test]
fn extract_uniprot_masked_low_memory() -> Result<()> {
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
            ">sp|P26709|1107L_ASFL5:128-138 RNELNNEEA 0",
            "RKELKKDEF%",
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
    for memory in &["", "-l", "-v"] {
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
    let mut args = vec!["locate".to_string(), filename.to_string()];

    if opts.absolute {
        args.push("-a".to_string());
    }

    if opts.low_memory {
        args.push("-l".to_string());
    }

    if let Some(max_query_len) = opts.max_query_len {
        args.push("-m".to_string());
        args.push(max_query_len.to_string());
    }

    for query in opts.queries {
        args.push(query.to_string());
    }

    let output = Command::cargo_bin(PRG)?.args(&args).output().expect("fail");
    assert!(output.status.success());

    let actual = String::from_utf8(output.stdout).expect("invalid UTF-8");
    let expected = fs::read_to_string(expected_file)?;
    assert_eq!(actual, expected);

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
            low_memory: false,
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
            low_memory: false,
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
            low_memory: false,
            max_query_len: None,
        },
        "../data/expected/uniprot-search1.out",
    )
}

// --------------------------------------------------
#[test]
fn locate_uniprot_masked() -> Result<()> {
    // cargo run -- lo data/expected/uniprot-masked.sufr RNEL DTPT GSGL
    // RNEL
    // sp|P0DO85|10H_STRNX 151
    // sp|P19084|11S3_HELAN 355
    // sp|P26709|1107L_ASFL5 128
    // sp|P32234|128UP_DROME 38
    // sp|Q6GZV8|017L_FRG3G 18
    // sp|Q6GZW6|009L_FRG3G 653
    // sp|Q9U408|14331_ECHGR 42
    // //
    // DTPT
    // sp|O55705|061R_IIV6 72
    // sp|Q196U6|114L_IIV3 123
    // sp|Q197D0|030L_IIV3 51
    // sp|Q197F8|002R_IIV3 368
    // sp|Q4U9M9|104K_THEAN 751
    // sp|Q6GZN2|093L_FRG3G 1
    // sp|Q6GZN6|089R_FRG3G 125
    // sp|Q6GZN8|087L_FRG3G 546
    // sp|Q6GZX3|002L_FRG3G 222
    // //
    // GSGL
    // sp|Q196Y8|072L_IIV3 29
    // sp|Q197B6|044L_IIV3 19
    // sp|Q197C3|037L_IIV3 40
    // sp|Q197D8|022L_IIV3 164
    // sp|Q67475|055R_FRG3G 119
    // sp|Q6GZR7|058R_FRG3G 150
    // sp|Q6GZT9|037R_FRG3G 119
    // sp|Q6GZU9|027R_FRG3G 496
    // sp|Q6GZV5|020R_FRG3G 115
    // sp|Q6GZV6|019R_FRG3G 816
    // sp|Q6GZW6|009L_FRG3G 394
    // sp|Q91G63|034R_IIV6 97
    // sp|Q9XHP0|11S2_SESIN 98
    // //

    locate(
        "../data/expected/uniprot-masked.sufr",
        LocateOptions {
            queries: vec!["RNEL".to_string(), "DTPT".to_string(), "GSGL".to_string()],
            absolute: false,
            low_memory: false,
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
            low_memory: false,
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
            low_memory: false,
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
            low_memory: true,
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
            low_memory: true,
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
        if let Some(caps) = regex.captures(&stdout) {
            assert_eq!(caps.get(1).unwrap().as_str().trim(), value);
        }
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
            ("Num Suffixes", "9"),
            ("Max query len", "0"),
            ("Seed mask", "None"),
            ("Num sequences", "1"),
            ("Sequence starts", "0"),
            ("Sequence names", "1"),
        ],
    )
}

// --------------------------------------------------
fn file_is_sorted(filename: &str, mask: Option<&str>) -> Result<()> {
    // Create the sufr file
    let sufr_file = NamedTempFile::new()?;
    let sufr_path = &sufr_file.path().to_string_lossy();
    let mut args = vec![
        "create".to_string(),
        "-o".to_string(),
        sufr_path.to_string(),
        filename.to_string(),
    ];

    let seed_mask: Option<SeedMask> = mask.map(SeedMask::new).transpose()?;
    if let Some(mask) = &seed_mask {
        args.push("-s".to_string());
        args.push(mask.to_string());
    }

    let output = Command::cargo_bin(PRG)?.args(&args).output()?;
    assert!(output.status.success());
    assert!(sufr_file.path().exists());

    // Use "list" to write the sorted suffixes to a temp file
    let list_file = NamedTempFile::new()?;
    let list_path = &list_file.path().to_string_lossy();
    let args = vec![
        "list".to_string(),
        "-o".to_string(),
        list_path.to_string(),
        sufr_path.to_string(),
    ];

    let output = Command::cargo_bin(PRG)?.args(&args).output()?;
    assert!(output.status.success());
    assert!(list_file.path().exists());

    // Read the lines of the sorted suffixes and ensure they are sorted
    let list_file = BufReader::new(File::open(list_file.path())?);
    let mut prev_line: Option<String> = None;
    for mut line in list_file.lines().map_while(Result::ok) {
        if let Some(mask) = &seed_mask {
            let chars: Vec<char> = line.chars().collect();
            let masked: String = mask
                .positions
                .iter()
                .filter_map(|&p| chars.get(p))
                .collect();
            line = masked;
        }

        if let Some(ref prev) = prev_line {
            assert!(*prev <= line);
        }

        prev_line = Some(line.clone());
    }

    Ok(())
}

// --------------------------------------------------
#[test]
fn long_dna_is_sorted() -> Result<()> {
    file_is_sorted(LONG, None)
}

// --------------------------------------------------
#[test]
fn masked_long_dna_is_sorted() -> Result<()> {
    file_is_sorted(LONG, Some("111010010100110111"))
}

// --------------------------------------------------
#[test]
fn uniprot_is_sorted() -> Result<()> {
    file_is_sorted(UNIPROT, None)
}

// --------------------------------------------------
#[test]
fn masked_uniprot_is_sorted() -> Result<()> {
    file_is_sorted(UNIPROT, Some("10111011"))
}
