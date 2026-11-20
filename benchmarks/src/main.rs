//! Benchmarks for `sufr`: each operation as a whole process, one row per
//! binary, drawn on a [`board`] as the measurements land.
//!
//! ```text
//! cargo run --release --                                     every benchmark
//! cargo run --release -- long-n -r 9                         one of them
//! cargo run --release -- -t 8 --cold                         8 threads, inputs read from disk
//! cargo run --release -- --compare-to before=/path/to/an/older/sufr
//! ```
//!
//! Every measurement is a whole process, because `wait4` is the one call
//! that reports wall, user and system time and peak memory together. So
//! process startup lands inside every row, which is why the fixtures are
//! generated, 10 MiB unless `--size` says otherwise, rather than the test
//! files under `data/`.
//!
//! The fixtures are generated once into `runs/` and kept; `--fresh` makes
//! them again. A measured run sends its stdout to `/dev/null` and its stderr
//! to `runs/stderr/`. Every measurement is also appended as one row to
//! `runs/results.tsv`, keyed by the run it belonged to, so two runs can be
//! joined later.
//!
//! With more than one `--bin`, the rows of each are interleaved and the order
//! flips each repeat, so a drift in the machine lands on both. Before the
//! measured runs, one unrecorded run of every row fingerprints its output,
//! and the run stops if two binaries disagree.

use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{SystemTime, UNIX_EPOCH};

use anyhow::{Context, bail};
use clap::Parser;

use bench::{BENCHES, Bench};
use board::{Board, Group};
use process::Timing;

mod bench;
mod board;
mod process;
mod synth;

// ---

#[derive(Parser)]
#[command(about = "sufr benchmarks")]
struct Cli {
    /// Benchmark to run, by name or prefix such as create; default all
    which: Option<String>,

    /// Threads for each sufr run
    #[arg(short, long, value_name = "N", default_value_t = 4)]
    threads: usize,

    /// A sufr binary as label=path; repeatable; default builds this tree
    #[arg(short, long = "bin", value_name = "LABEL=PATH")]
    bins: Vec<String>,

    /// Another sufr binary, as a second row beside this tree's
    #[arg(long, value_name = "LABEL=PATH", conflicts_with = "bins")]
    compare_to: Option<String>,

    /// Runs per row
    #[arg(short, long, value_name = "N", default_value_t = 5)]
    repeat: usize,

    /// Fixture size in MiB
    #[arg(short, long, value_name = "MIB", default_value_t = 10)]
    size: usize,

    /// Evict the input from the page cache before each run
    #[arg(short, long)]
    cold: bool,

    /// Directory for fixtures, results and stderr
    #[arg(short, long, value_name = "DIR")]
    out: Option<PathBuf>,

    /// Regenerate the fixtures
    #[arg(short, long)]
    fresh: bool,
}

/// One sufr binary on the board.
struct Bin {
    label: String,
    path: PathBuf,
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    let chosen: Vec<&Bench> = match &cli.which {
        None => BENCHES.iter().collect(),
        Some(which) => {
            let found: Vec<&Bench> = BENCHES
                .iter()
                .filter(|b| b.name == which || b.name.starts_with(&format!("{which}-")))
                .collect();
            if found.is_empty() {
                let known: Vec<&str> = BENCHES.iter().map(|b| b.name).collect();
                bail!(
                    "no benchmark called {which:?}; there is {}",
                    known.join(", ")
                );
            }
            found
        }
    };

    let out = cli
        .out
        .clone()
        .unwrap_or_else(|| Path::new(env!("CARGO_MANIFEST_DIR")).join("runs"));
    fs::create_dir_all(out.join("stderr"))
        .with_context(|| format!("creating {}", out.display()))?;

    let bins = bins(&cli.bins, cli.compare_to.as_deref())?;
    for bin in &bins {
        println!("{} = {}", bin.label, bin.path.display());
    }

    let fixture = Fixture::at(&out, cli.size);
    fixture.ensure(cli.fresh, cli.size << 20)?;

    let header = format!(
        "{} MiB fixtures · {} threads · median of {} · {}",
        cli.size,
        cli.threads,
        cli.repeat,
        if cli.cold { "cold" } else { "warm" },
    );

    let plan: Vec<Group> = chosen
        .iter()
        .map(|bench| Group {
            op: bench.name,
            why: bench.same_work,
            rows: bins
                .iter()
                .enumerate()
                .map(|(n, bin)| (bin.label.clone(), n))
                .collect(),
        })
        .collect();

    let mut results = Results::open(&out)?;
    let run = Run {
        bins: &bins,
        fixture: &fixture,
        out: &out,
        threads: cli.threads,
        repeat: cli.repeat,
        cold: cli.cold,
    };

    let board = Board::new(&plan, header, cli.repeat, bins[0].label.clone());
    let painter = board.start();

    let ran = chosen
        .iter()
        .try_for_each(|bench| run.one(bench, &board, &mut results));

    // the thread lets go of stdout before anything else
    // writes to it, including the error that stopped the run
    painter.stop();
    ran?;

    board.finish();

    Ok(())
}

/// The binaries from `--bin label=path`, or this tree built and labelled
/// `head`, with the `--compare-to` binary after it.
fn bins(given: &[String], compare_to: Option<&str>) -> anyhow::Result<Vec<Bin>> {
    if given.is_empty() {
        let mut bins = vec![build()?];
        if let Some(spec) = compare_to {
            bins.push(parse_bin(spec)?);
        }
        return Ok(bins);
    }

    given.iter().map(|spec| parse_bin(spec)).collect()
}

/// `label=path`, or a path labelled by its file name.
fn parse_bin(spec: &str) -> anyhow::Result<Bin> {
    let (label, path) = match spec.split_once('=') {
        Some((label, path)) => (label.to_string(), PathBuf::from(path)),
        None => {
            let path = PathBuf::from(spec);
            let label = path
                .file_name()
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_default();
            (label, path)
        }
    };
    if !path.is_file() {
        bail!("{spec}: {} is not a file", path.display());
    }
    Ok(Bin { label, path })
}

/// `cargo build --release` in the tree above, and the binary it leaves.
fn build() -> anyhow::Result<Bin> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    println!("building sufr ...");

    let status = Command::new(option_env!("CARGO").unwrap_or("cargo"))
        .current_dir(&root)
        .args(["build", "--release", "--bin", "sufr"])
        .status()
        .context("running cargo")?;
    if !status.success() {
        bail!("cargo build failed");
    }

    Ok(Bin {
        label: "head".to_string(),
        path: root.join("target/release/sufr"),
    })
}

// --- the run

struct Run<'a> {
    bins: &'a [Bin],
    fixture: &'a Fixture,
    out: &'a Path,
    threads: usize,
    repeat: usize,
    cold: bool,
}

impl Run<'_> {
    fn one(
        &self,
        bench: &Bench,
        board: &Board,
        results: &mut Results,
    ) -> anyhow::Result<()> {
        let scratch = self.out.join("scratch.sufr");
        let dry_out = self.out.join("dry.out");

        // one unrecorded run of each row, so a measured one
        // is never the first to load the binary and never
        // the first to fail. its output is what the binaries
        // are checked against each other on
        let mut seen: Option<(&str, u64)> = None;
        for bin in self.bins {
            remove(&scratch)?;
            let mut cmd = self.command(bin, bench, &scratch);
            board.running(
                bench.name,
                &bin.label,
                format!("dry {}", bin.label),
                line(&cmd),
            );

            cmd.stdout(File::create(&dry_out)?);
            process::measure(&mut cmd)?;
            wrote(bench, bin, &scratch)?;
            let print = process::fingerprint(&[&dry_out, &scratch])?;

            match seen {
                None => seen = Some((&bin.label, print)),
                Some((first, p)) if p != print => bail!(
                    "{}: {first} and {} write different output",
                    bench.name,
                    bin.label
                ),
                Some(_) => {}
            }
        }
        remove(&dry_out)?;

        for rep in 1..=self.repeat {
            // the order flips each repeat, so a machine
            // that slows down over the run slows both
            let order: Vec<&Bin> = if rep % 2 == 1 {
                self.bins.iter().collect()
            } else {
                self.bins.iter().rev().collect()
            };

            for bin in order {
                remove(&scratch)?;
                if self.cold {
                    process::evict(&[bench.reads(self.fixture)])?;
                }

                let mut cmd = self.command(bin, bench, &scratch);
                board.running(
                    bench.name,
                    &bin.label,
                    format!("{} #{rep}", bin.label),
                    line(&cmd),
                );

                cmd.stdout(Stdio::null());
                let timing = process::measure(&mut cmd)?;
                wrote(bench, bin, &scratch)?;

                board.record(bench.name, &bin.label, timing);
                results.record(
                    &bin.label,
                    bench.name,
                    rep,
                    self.threads,
                    self.cold,
                    timing,
                )?;
            }
        }

        remove(&scratch)
    }

    fn command(&self, bin: &Bin, bench: &Bench, scratch: &Path) -> Command {
        let stderr = self
            .out
            .join("stderr")
            .join(format!("{}-{}.txt", bench.name, bin.label));

        // sufr takes -t after the subcommand
        let argv = bench.argv(self.fixture, scratch);
        let mut cmd = Command::new(&bin.path);
        cmd.arg(&argv[0])
            .arg("-t")
            .arg(self.threads.to_string())
            .args(&argv[1..])
            .stderr(File::create(stderr).expect("runs/stderr exists"));
        cmd
    }
}

/// The command as one line, for the status line under the board.
fn line(cmd: &Command) -> String {
    std::iter::once(cmd.get_program())
        .chain(cmd.get_args())
        .map(|a| a.to_string_lossy())
        .collect::<Vec<_>>()
        .join(" ")
}

/// An error unless `bin` left an index at `scratch`.
fn wrote(bench: &Bench, bin: &Bin, scratch: &Path) -> anyhow::Result<()> {
    if !scratch.is_file() {
        bail!(
            "{}: {} exited 0 but wrote no index at {}",
            bench.name,
            bin.label,
            scratch.display()
        );
    }
    Ok(())
}

fn remove(path: &Path) -> anyhow::Result<()> {
    match fs::remove_file(path) {
        Ok(()) => Ok(()),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(e) => Err(e).with_context(|| format!("removing {}", path.display())),
    }
}

// --- the fixtures

/// A fixture generator: the path to write and the size in bases.
type Writer = fn(&Path, usize) -> std::io::Result<()>;

/// The generated inputs of one size, kept between runs.
pub(crate) struct Fixture {
    pub(crate) dna: PathBuf,
    pub(crate) long_n: PathBuf,
    pub(crate) long_repeats: PathBuf,
    pub(crate) protein: PathBuf,
}

impl Fixture {
    fn at(out: &Path, size: usize) -> Fixture {
        Fixture {
            dna: out.join(format!("dna-{size}.fa")),
            long_n: out.join(format!("long-n-{size}.fa")),
            long_repeats: out.join(format!("long-repeats-{size}.fa")),
            protein: out.join(format!("protein-{size}.fa")),
        }
    }

    /// Generate whatever is missing, or everything when `fresh`.
    fn ensure(&self, fresh: bool, bases: usize) -> anyhow::Result<()> {
        let writers: [(&Path, Writer); 4] = [
            (&self.dna, synth::dna),
            (&self.long_n, synth::long_n),
            (&self.long_repeats, synth::long_repeats),
            (&self.protein, synth::protein),
        ];

        for (path, write) in writers {
            if fresh {
                remove(path)?;
            }
            if !path.is_file() {
                println!("generating {} ...", path.display());
                write(path, bases)?;
            }
        }

        Ok(())
    }
}

// --- the results

/// `runs/results.tsv`, one row per measurement, appended to across runs.
struct Results {
    file: File,

    /// Seconds since the epoch when this run started, which is what joins
    /// its rows.
    run: u64,
}

impl Results {
    fn open(out: &Path) -> anyhow::Result<Results> {
        let path = out.join("results.tsv");
        let new = !path.is_file();

        let mut file = File::options()
            .append(true)
            .create(true)
            .open(&path)
            .with_context(|| format!("opening {}", path.display()))?;
        if new {
            writeln!(
                file,
                "run\tbin\tbench\trep\tthreads\tcold\twall_s\tuser_s\tsys_s\tmax_rss_kb"
            )?;
        }

        Ok(Results {
            file,
            run: SystemTime::now().duration_since(UNIX_EPOCH)?.as_secs(),
        })
    }

    fn record(
        &mut self,
        bin: &str,
        bench: &str,
        rep: usize,
        threads: usize,
        cold: bool,
        t: Timing,
    ) -> anyhow::Result<()> {
        writeln!(
            self.file,
            "{}\t{bin}\t{bench}\t{rep}\t{threads}\t{}\t{:.4}\t{:.4}\t{:.4}\t{}",
            self.run,
            if cold { "cold" } else { "warm" },
            t.wall_s,
            t.user_s,
            t.sys_s,
            t.max_rss_kb
        )?;
        Ok(())
    }
}
