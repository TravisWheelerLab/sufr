use anyhow::Result;
use chrono::{DateTime, Local};
use clap::{builder::PossibleValue, Parser, ValueEnum};
use format_num::NumberFormat;
use libsufr::{
    read_sequence_file, read_suffix_length, FromUsize, Int, SufrBuilder,
    SufrBuilderArgs, SufrFile,
};
use log::info;
use regex::Regex;
use std::{
    cmp::min,
    ffi::OsStr,
    fmt::Debug,
    fs::{self, File},
    io::{self, Write},
    path::{Path, PathBuf},
    time::Instant,
};
use tabled::Table;

// --------------------------------------------------
#[derive(Parser, Debug)]
#[command(arg_required_else_help = true, version, about)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Option<Command>,

    /// Log level
    #[arg(short, long)]
    pub log: Option<LogLevel>,

    /// Log file
    #[arg(long)]
    pub log_file: Option<String>,
}

#[derive(Parser, Debug)]
pub enum Command {
    /// Create sufr file
    Create(CreateArgs),

    /// Check sufr file for correctness
    Check(CheckArgs),

    /// Extract suffixes from a sufr file
    Extract(ExtractArgs),

    /// Locate for sequences in a sufr file
    Locate(LocateArgs),

    /// Summarize sufr file
    Summarize(SummarizeArgs),
}

#[derive(Debug, Parser)]
#[command(about, alias = "ch")]
pub struct CheckArgs {
    /// Sufr file
    #[arg(value_name = "SUFR")]
    pub filename: String,

    /// List errors
    #[arg(short, long)]
    pub verbose: bool,
}

#[derive(Debug, Parser)]
#[command(about, alias = "cr")]
pub struct CreateArgs {
    /// Input file
    #[arg(value_name = "INPUT")]
    pub input: String,

    /// Subproblem count
    #[arg(short, long, value_name = "NUM_PARTS", default_value = "16")]
    pub num_partitions: usize,

    /// Max context
    #[arg(short, long, value_name = "CONTEXT")]
    pub max_query_len: Option<usize>,

    /// Number of threads
    #[arg(short, long, value_name = "THREADS")]
    pub threads: Option<usize>,

    /// Output file
    #[arg(short, long, value_name = "OUTPUT")]
    pub output: Option<String>,

    /// Input is DNA
    #[arg(short('d'), long("dna"))]
    pub is_dna: bool,

    /// Allow suffixes starting with ambiguity codes
    #[arg(short, long)]
    pub allow_ambiguity: bool,

    /// Ignore suffixes in soft-mask/lowercase regions
    #[arg(short, long)]
    pub ignore_softmask: bool,
}

#[derive(Debug, Parser)]
#[command(about, alias = "ex")]
pub struct ExtractArgs {
    /// Maximum length of sequence
    #[arg(short, long, value_name = "MAX")]
    pub max_len: Option<usize>,

    /// Output
    #[arg(short, long)]
    pub output: Option<String>,

    /// Sufr file
    #[arg(value_name = "SUFR")]
    pub filename: String,

    /// Suffixes to extract
    #[arg(
        value_name = "SUFFIX", 
        value_parser = clap::value_parser!(u64).range(0..),
        num_args(1..),
    )]
    pub suffixes: Vec<u64>,
}

#[derive(Debug, Parser)]
#[command(about, alias = "lo")]
pub struct LocateArgs {
    /// Maximum query length
    #[arg(short, long, value_name = "LEN")]
    pub max_query_len: Option<usize>,

    /// Output
    #[arg(short, long, value_name = "OUT")]
    pub output: Option<String>,

    /// Sufr file
    #[arg(value_name = "SUFR")]
    pub file: String,

    /// Query
    #[arg(value_name = "QUERY", required = true)]
    pub query: Vec<String>,
}

#[derive(Debug, Parser)]
#[command(about, alias = "su")]
pub struct SummarizeArgs {
    /// Sufr file
    #[arg(value_name = "SUFR")]
    pub file: String,
}

#[derive(Debug, Clone)]
pub enum LogLevel {
    Info,
    Debug,
}

impl ValueEnum for LogLevel {
    fn value_variants<'a>() -> &'a [Self] {
        &[LogLevel::Info, LogLevel::Debug]
    }

    fn to_possible_value<'a>(&self) -> Option<PossibleValue> {
        Some(match self {
            LogLevel::Info => PossibleValue::new("info"),
            LogLevel::Debug => PossibleValue::new("debug"),
        })
    }
}

// --------------------------------------------------
pub fn check(args: &CheckArgs) -> Result<()> {
    let len = read_suffix_length(&args.filename)? as u64;
    if len < u32::MAX as u64 {
        let sufr_file: SufrFile<u32> = SufrFile::read(&args.filename)?;
        _check(sufr_file, args)
    } else {
        let sufr_file: SufrFile<u64> = SufrFile::read(&args.filename)?;
        _check(sufr_file, args)
    }
}

// --------------------------------------------------
pub fn _check<T>(mut sufr_file: SufrFile<T>, _args: &CheckArgs) -> Result<()>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + Debug,
{
    let now = Instant::now();
    let errors = sufr_file.check_suffix_array()?;
    let num_errors = errors.len();
    let num_fmt = NumberFormat::new();
    println!(
        "Checked {} suffix{}, found {} error{} in suffix array.",
        num_fmt.format(",.0", sufr_file.num_suffixes.to_usize() as f64),
        if sufr_file.num_suffixes.to_usize() == 1 {
            ""
        } else {
            "es"
        },
        num_fmt.format(",.0", num_errors as f64),
        if num_errors == 1 { "" } else { "s" },
    );

    //let errors = sa.check_order();
    //let num_errors = errors.len();
    //
    //if args.verbose {
    //    for (i, pos) in errors.iter().enumerate() {
    //        println!(
    //            "{:3}: pos {pos:5} {}",
    //            i + 1,
    //            sa.string_at(pos.to_usize())
    //        );
    //    }
    //}
    //
    //let lcp_errors = sa.check_lcp();
    //let num_errors = lcp_errors.len();
    //println!(
    //    "Found {num_errors} error{} in LCP",
    //    if num_errors == 1 { "" } else { "s" },
    //);
    //
    //if args.verbose {
    //    for (i, pos) in lcp_errors.iter().enumerate() {
    //        println!("{:3}: pos {pos:5}", i + 1,);
    //    }
    //}
    //
    println!("Finished checking in {:?}.", now.elapsed());
    Ok(())
}

// --------------------------------------------------
pub fn create(args: &CreateArgs) -> Result<()> {
    let num_threads = args.threads.unwrap_or(num_cpus::get());
    info!(
        "Using {num_threads} thread{}",
        if num_threads == 1 { "" } else { "s" }
    );
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .unwrap();

    // Read sequence input
    let now = Instant::now();
    let sequence_delimiter = if args.is_dna { b'N' } else { b'X' };
    let seq_data = read_sequence_file(&args.input, sequence_delimiter)?;
    let text_len = seq_data.seq.len();
    let num_fmt = NumberFormat::new();
    info!(
        "Read input of len {} in {:?}",
        num_fmt.format(",.0", text_len as f64),
        now.elapsed()
    );
    let builder_args = SufrBuilderArgs {
        text: seq_data.seq,
        max_query_len: args.max_query_len,
        is_dna: args.is_dna,
        allow_ambiguity: args.allow_ambiguity,
        ignore_softmask: args.ignore_softmask,
        sequence_starts: seq_data.start_positions.into_iter().collect(),
        headers: seq_data.headers,
        num_partitions: args.num_partitions,
        sequence_delimiter,
    };

    if (text_len as u64) < u32::MAX as u64 {
        let sa: SufrBuilder<u32> = SufrBuilder::new(builder_args)?;
        _create(sa, args, now)
    } else {
        let sa: SufrBuilder<u64> = SufrBuilder::new(builder_args)?;
        _create(sa, args, now)
    }
}

// --------------------------------------------------
// Helper for "create" that checks/writes SA/LCP
pub fn _create<T>(sa: SufrBuilder<T>, args: &CreateArgs, timer: Instant) -> Result<()>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + Debug,
{
    let outfile = &args.output.clone().unwrap_or(format!(
        "{}.sufr",
        PathBuf::from(&args.input)
            .file_stem()
            .unwrap_or(OsStr::new("out"))
            .to_string_lossy()
    ));
    let now = Instant::now();
    let bytes_written = sa.write(outfile)?;
    let num_fmt = NumberFormat::new();
    info!(
        "Wrote {} byte{} to '{outfile}' in {:?}",
        num_fmt.format(",.0", bytes_written as f64),
        if bytes_written == 1 { "" } else { "s" },
        now.elapsed()
    );
    info!("Total time: {:?}", timer.elapsed());

    Ok(())
}

// --------------------------------------------------
pub fn extract(args: &ExtractArgs) -> Result<()> {
    let len = read_suffix_length(&args.filename)? as u64;
    if len < u32::MAX as u64 {
        let now = Instant::now();
        let sufr_file: SufrFile<u32> = SufrFile::read(&args.filename)?;
        info!("Read sufr file in {:?}", now.elapsed());
        _extract(sufr_file, args)
    } else {
        let now = Instant::now();
        let sufr_file: SufrFile<u64> = SufrFile::read(&args.filename)?;
        info!("Read sufr file in {:?}", now.elapsed());
        _extract(sufr_file, args)
    }
}

// --------------------------------------------------
pub fn _extract<T>(sufr_file: SufrFile<T>, args: &ExtractArgs) -> Result<()>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + Debug,
{
    let now = Instant::now();

    let mut output: Box<dyn Write> = match &args.output {
        Some(out_name) => Box::new(File::create(out_name)?),
        _ => Box::new(io::stdout()),
    };

    let text_len = sufr_file.text_len.to_usize();
    for suffix in args.suffixes.iter().map(|&v| v as usize) {
        let end = min(args.max_len.map_or(text_len, |len| suffix + len), text_len);
        let seq = String::from_utf8(sufr_file.text[suffix..end].to_vec())?;
        writeln!(output, "{seq}")?;
    }
    info!("Extracted suffixes in {:?}", now.elapsed());

    Ok(())
}

// --------------------------------------------------
pub fn locate(args: &LocateArgs) -> Result<()> {
    let len = read_suffix_length(&args.file)? as u64;
    if len < u32::MAX as u64 {
        let now = Instant::now();
        let sa: SufrFile<u32> = SufrFile::read(&args.file)?;
        info!("Read sufr file in {:?}", now.elapsed());
        _locate(sa, args)
    } else {
        let now = Instant::now();
        let sa: SufrFile<u64> = SufrFile::read(&args.file)?;
        info!("Read sufr file in {:?}", now.elapsed());
        _locate(sa, args)
    }
}

// --------------------------------------------------
pub fn _locate<T>(mut sufr_file: SufrFile<T>, args: &LocateArgs) -> Result<()>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + Debug,
{
    let mut output: Box<dyn Write> = match &args.output {
        Some(out_name) => Box::new(File::create(out_name)?),
        _ => Box::new(io::stdout()),
    };

    let whitespace = Regex::new(r"\s+").unwrap();
    let mut queries = vec![];
    for query in &args.query {
        if Path::new(&query).exists() {
            let contents = fs::read_to_string(query)?;
            let mut vals: Vec<String> = whitespace
                .split(&contents)
                .filter(|v| !v.is_empty())
                .map(|v| v.to_string())
                .collect();
            queries.append(&mut vals);
        } else {
            queries.push(query.to_string());
        }
    }

    let now = Instant::now();
    for res in sufr_file.locate(&queries, args.max_query_len)? {
        let loc = if res.suffixes.is_empty() {
            "not found".to_string()
        } else {
            res.suffixes.iter().map(|s| s.to_string()).collect::<Vec<_>>().join(" ")
        };
        writeln!(output, "{}: {loc}", res.query)?;
    }

    info!(
        "Locate of {} finished in {:?}",
        queries.len(),
        now.elapsed()
    );
    Ok(())
}

// --------------------------------------------------
pub fn summarize(args: &SummarizeArgs) -> Result<()> {
    let len = read_suffix_length(&args.file)? as u64;
    if len < u32::MAX as u64 {
        let now = Instant::now();
        let sa: SufrFile<u32> = SufrFile::read(&args.file)?;
        info!("Read sufr file in {:?}", now.elapsed());
        _summarize(sa, args)
    } else {
        let now = Instant::now();
        let sa: SufrFile<u64> = SufrFile::read(&args.file)?;
        info!("Read sufr file in {:?}", now.elapsed());
        _summarize(sa, args)
    }
}

// --------------------------------------------------
pub fn _summarize<T>(sufr_file: SufrFile<T>, _args: &SummarizeArgs) -> Result<()>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + Debug,
{
    let num_fmt = NumberFormat::new();
    let mut rows = vec![vec!["Filename".to_string(), sufr_file.filename.to_string()]];
    let metadata = fs::metadata(&sufr_file.filename)?;
    let modified: DateTime<Local> = DateTime::from(metadata.modified()?);
    rows.push(vec![
        "Modified".to_string(),
        modified.format("%Y-%m-%d %H:%M").to_string(),
    ]);
    rows.push(vec![
        "File Size".to_string(),
        format!(
            "{} bytes",
            num_fmt.format(",.0", metadata.len().to_usize() as f64)
        ),
    ]);
    rows.push(vec![
        "File Version".to_string(),
        sufr_file.version.to_string(),
    ]);
    rows.push(vec!["DNA".to_string(), sufr_file.is_dna.to_string()]);
    rows.push(vec![
        "Text Length".to_string(),
        num_fmt.format(",.0", sufr_file.text_len.to_usize() as f64),
    ]);
    rows.push(vec![
        "Num Suffixes".to_string(),
        num_fmt.format(",.0", sufr_file.num_suffixes.to_usize() as f64),
    ]);
    rows.push(vec![
        "Max query len".to_string(),
        num_fmt.format(",.0", sufr_file.max_query_len.to_usize() as f64),
    ]);
    rows.push(vec![
        "Num sequences".to_string(),
        num_fmt.format(",.0", sufr_file.num_sequences.to_usize() as f64),
    ]);
    let seq_starts = sufr_file
        .sequence_starts
        .into_iter()
        .map(|v| v.to_string())
        .collect::<Vec<_>>()
        .join(", ");
    rows.push(vec![
        "Sequence starts".to_string(),
        textwrap::wrap(&seq_starts, 40).join("\n"),
    ]);
    let seq_headers = sufr_file.headers.into_iter().collect::<Vec<_>>().join(", ");
    rows.push(vec![
        "Headers".to_string(),
        textwrap::wrap(&seq_headers, 40).join("\n"),
    ]);

    let table = Table::from_iter(rows);
    println!("{table}");

    Ok(())
}
