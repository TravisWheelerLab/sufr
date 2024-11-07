use anyhow::Result;
use assert_cmd::Command;
use std::fs;
use tempfile::NamedTempFile;

const PRG: &str = "sufr";
const SEQ1: &str = "tests/inputs/1.fa";
const SEQ2: &str = "tests/inputs/2.fa";
const SEQ3: &str = "tests/inputs/3.fa";

struct CreateArgs {
    allow_ambiguity: bool,
    ignore_softmask: bool,
}

// --------------------------------------------------
fn create(input_file: &str, expected_file: &str, create_args: CreateArgs) -> Result<()> {
    let outfile = NamedTempFile::new()?;
    let outpath = &outfile.path().to_str().unwrap();
    let mut args = vec!["create", "--dna", "-o", outpath, input_file];
    if create_args.allow_ambiguity {
        args.push("--allow-ambiguity");
    }
    if create_args.ignore_softmask {
        args.push("--ignore-softmask");
    }
    let output = Command::cargo_bin(PRG)?.args(args).output().expect("fail");

    assert!(output.status.success());
    assert!(outfile.path().exists());

    let actual = fs::read(outfile.path())?;
    let expected = fs::read(expected_file)?;

    assert_eq!(actual, expected);

    Ok(())
}

// --------------------------------------------------
fn check(filename: &str) -> Result<()> {
    let output = Command::cargo_bin(PRG)?
        .args(&["check", filename])
        .output()
        .expect("fail");

    assert!(output.status.success());

    let stdout = String::from_utf8(output.stdout).expect("invalid UTF-8");

    assert!(stdout.contains("found 0 errors"));

    Ok(())
}

// --------------------------------------------------
#[test]
fn create_empty_dies() -> Result<()> {
    Command::cargo_bin(PRG)?
        .args(["create", "tests/expected/empty.sa"])
        .assert()
        .failure();
    Ok(())
}

// --------------------------------------------------
#[test]
fn create_seq1() -> Result<()> {
    let args = CreateArgs { allow_ambiguity: false, ignore_softmask: false };
    create(SEQ1, "tests/expected/1.sufr", args)
}

// --------------------------------------------------
#[test]
fn create_seq1_allow_ambiguity() -> Result<()> {
    let args = CreateArgs { allow_ambiguity: true, ignore_softmask: false };
    create(SEQ1, "tests/expected/1n.sufr", args)
}

// --------------------------------------------------
#[test]
fn create_seq2() -> Result<()> {
    let args = CreateArgs { allow_ambiguity: false, ignore_softmask: false };
    create(SEQ2, "tests/expected/2.sufr", args)
}

// --------------------------------------------------
#[test]
fn create_seq2_allow_ambiguity() -> Result<()> {
    let args = CreateArgs { allow_ambiguity: true, ignore_softmask: false };
    create(SEQ2, "tests/expected/2n.sufr", args)
}

// --------------------------------------------------
#[test]
fn create_seq2_ignore_softmask() -> Result<()> {
    let args = CreateArgs { allow_ambiguity: false, ignore_softmask: true };
    create(SEQ2, "tests/expected/2s.sufr", args)
}

// --------------------------------------------------
#[test]
fn create_seq2_allow_ambiguity_ignore_softmask() -> Result<()> {
    let args = CreateArgs { allow_ambiguity: true, ignore_softmask: true };
    create(SEQ2, "tests/expected/2ns.sufr", args)
}

// --------------------------------------------------
#[test]
fn create_seq3() -> Result<()> {
    let args = CreateArgs { allow_ambiguity: false, ignore_softmask: false };
    create(SEQ3, "tests/expected/3.sufr", args)
}

// --------------------------------------------------
#[test]
fn check_seq1() -> Result<()> {
    check("tests/expected/1.sufr")
}

// --------------------------------------------------
//#[test]
//fn read_input1_1_10() -> Result<()> {
//    read(SEQ1, SEQ1_SA, "1-10", "tests/expected/seq1.read.1-10.out")
//}

// --------------------------------------------------
//#[test]
//fn read_input1_1_5() -> Result<()> {
//    read(SEQ1, SEQ1_SA, "1-5", "tests/expected/seq1.read.1-5.out")
//}

// --------------------------------------------------
//#[test]
//fn read_input1_5_10() -> Result<()> {
//    read(SEQ1, SEQ1_SA, "5-10", "tests/expected/seq1.read.5-10.out")
//}

// --------------------------------------------------
//#[test]
//fn check_input1() -> Result<()> {
//    check(SEQ1, SEQ1_SA)
//}
