//! What is measured: one row group per fixture, with a row for each binary.

use std::ffi::OsString;
use std::path::{Path, PathBuf};

use crate::Fixture;

/// One row group.
pub(crate) struct Bench {
    /// The group label, and what `WHICH` on the command line matches.
    pub(crate) name: &'static str,

    /// What every row of the group does, printed under the board. A reader
    /// who cannot check the claim cannot read the numbers.
    pub(crate) same_work: &'static str,

    pub(crate) flags: &'static [&'static str],
    pub(crate) input: Input,
}

pub(crate) enum Input {
    Dna,
    LongN,
    LongRepeats,
    Protein,
}

pub(crate) const BENCHES: &[Bench] = &[
    Bench {
        name: "normal",
        same_work: "DNA with diverged copies of 50 repeat families; the control",
        flags: &["--dna", "--write-lcp"],
        input: Input::Dna,
    },
    Bench {
        name: "long-n",
        same_work: "the same DNA with runs of N inserted, 1%, 3% and 10% of the size, N indexed over the byte alphabet",
        flags: &["--write-lcp"],
        input: Input::LongN,
    },
    Bench {
        name: "long-repeats",
        same_work: "a tenth of the DNA ten times over, one random base between copies",
        flags: &["--dna", "--write-lcp"],
        input: Input::LongRepeats,
    },
    Bench {
        name: "protein",
        same_work: "protein with log-normal record lengths, over the byte alphabet",
        flags: &[],
        input: Input::Protein,
    },
    Bench {
        name: "mask-dna",
        same_work: "the normal DNA under the seed mask 11011011011011011, 12 care positions of 17",
        flags: &["--dna", "-s", "11011011011011011"],
        input: Input::Dna,
    },
    Bench {
        name: "mask-dna-heavy",
        same_work: "the normal DNA under the seed mask 11011011011011011011011011011, 20 care positions of 29",
        flags: &["--dna", "-s", "11011011011011011011011011011"],
        input: Input::Dna,
    },
    Bench {
        name: "mask-dna-wide",
        same_work: "the normal DNA under a seed mask with a 1 at every fourth position, 21 care positions of 81",
        flags: &[
            "--dna",
            "-s",
            "100010001000100010001000100010001000100010001000100010001000100010001000100010001",
        ],
        input: Input::Dna,
    },
    Bench {
        name: "mask-protein",
        same_work: "the protein under the seed mask 11011011011011011, 12 care positions of 17",
        flags: &["-s", "11011011011011011"],
        input: Input::Protein,
    },
];

impl Bench {
    /// The arguments after `sufr -t N`, writing to `scratch`.
    pub(crate) fn argv(&self, fixture: &Fixture, scratch: &Path) -> Vec<OsString> {
        let mut argv: Vec<OsString> = vec!["create".into()];
        argv.extend(self.flags.iter().map(OsString::from));
        argv.extend(["-o".into(), scratch.into(), self.reads(fixture).into()]);
        argv
    }

    /// What a run reads, which is what `--cold` evicts.
    pub(crate) fn reads(&self, fixture: &Fixture) -> PathBuf {
        match self.input {
            Input::Dna => fixture.dna.clone(),
            Input::LongN => fixture.long_n.clone(),
            Input::LongRepeats => fixture.long_repeats.clone(),
            Input::Protein => fixture.protein.clone(),
        }
    }
}
