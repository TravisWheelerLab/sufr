use anyhow::{anyhow, bail, Result};
use chrono::{DateTime, Local};
use clap::{builder::PossibleValue, Parser, ValueEnum};
use format_num::NumberFormat;
use libsufr::{
    sufr_builder::{SufrBuilder, SufrBuilderArgs},
    sufr_file::SufrFile,
    types::{FromUsize, Int, SearchOptions},
    util::{read_sequence_file, read_text_length},
};
use log::info;
use regex::Regex;
use std::{
    cmp::min,
    ffi::OsStr,
    fmt::Debug,
    fs::{self, File},
    io::{self, Write},
    ops::Range,
    path::{Path, PathBuf},
    time::Instant,
};
use tabled::Table;

// --------------------------------------------------
type PositionList = Vec<Range<usize>>;

// --------------------------------------------------
#[derive(Parser, Debug)]
#[command(arg_required_else_help = true, version, about)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Option<Command>,

    /// Number of threads
    #[arg(short, long, value_name = "THREADS")]
    pub threads: Option<usize>,

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

    /// List the suffix array from a sufr file
    List(ListArgs),

    /// Count occurrences of sequences in a sufr file
    Count(CountArgs),

    /// Locate sequences in a sufr file
    Locate(LocateArgs),

    /// Summarize sufr file
    Summarize(SummarizeArgs),
}

#[derive(Debug, Parser)]
#[command(about, alias = "ch")]
pub struct CheckArgs {
    /// Sufr file
    #[arg(value_name = "SUFR")]
    pub file: String,

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

    /// Character to separate sequences
    #[arg(short('D'), long, default_value = "%", value_name = "DELIM")]
    pub sequence_delimiter: char,
}

#[derive(Debug, Parser)]
#[command(about, alias = "ex")]
pub struct ExtractArgs {
    /// Prefix length
    #[arg(short, long, value_name = "PREFIX_LEN")]
    pub prefix_len: Option<usize>,

    /// Suffix length
    #[arg(short, long, value_name = "SUFFIX_LEN")]
    pub suffix_len: Option<usize>,

    /// Show LCP
    #[arg(short, long)]
    pub lcp: bool,

    /// Output
    #[arg(short, long, value_name = "OUT")]
    pub output: Option<String>,

    /// Sufr file
    #[arg(value_name = "SUFR")]
    pub file: String,

    /// Suffixes to extract
    #[arg(value_name = "SUFFIX", num_args(1..))]
    pub suffixes: Vec<String>,
}

#[derive(Debug, Parser)]
#[command(about, alias = "ls")]
pub struct ListArgs {
    /// Length of suffixes to show
    #[arg(short, long, value_name = "LEN")]
    pub len: Option<usize>,

    /// Output
    #[arg(short, long, value_name = "OUT")]
    pub output: Option<String>,

    /// Sufr file
    #[arg(value_name = "SUFR")]
    pub file: String,

    /// Ranks of suffixes to show
    #[arg(value_name = "POS")]
    pub ranks: Vec<String>,
}

#[derive(Debug, Parser)]
#[command(about, alias = "co")]
pub struct CountArgs {
    /// Maximum query length
    #[arg(short, long, value_name = "LEN")]
    pub max_query_len: Option<usize>,

    /// Output
    #[arg(short, long, value_name = "OUT")]
    pub output: Option<String>,

    /// Low-memory
    #[arg(short, long)]
    pub low_memory: bool,

    /// Sufr file
    #[arg(value_name = "SUFR")]
    pub file: String,

    /// Query
    #[arg(value_name = "QUERY", required = true)]
    pub query: Vec<String>,
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

    /// Low-memory
    #[arg(short, long)]
    pub low_memory: bool,

    /// Show absolute position in text
    #[arg(short, long)]
    pub abs: bool,

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

// --------------------------------------------------
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
    let text_len = read_text_length(&args.file)? as u64;
    if text_len < u32::MAX as u64 {
        let sufr_file: SufrFile<u32> = SufrFile::read(&args.file)?;
        _check(sufr_file, args)
    } else {
        let sufr_file: SufrFile<u64> = SufrFile::read(&args.file)?;
        _check(sufr_file, args)
    }
}

// --------------------------------------------------
fn _check<T>(mut sufr_file: SufrFile<T>, args: &CheckArgs) -> Result<()>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + Debug,
{
    let now = Instant::now();
    let errors = sufr_file.check()?;
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

    if args.verbose {
        for err in errors {
            println!("{err}");
        }
    }
    println!("Finished checking in {:?}.", now.elapsed());
    Ok(())
}

// --------------------------------------------------
pub fn count(args: &CountArgs) -> Result<()> {
    let text_len = read_text_length(&args.file)? as u64;
    if text_len < u32::MAX as u64 {
        let sufr_file: SufrFile<u32> = SufrFile::read(&args.file)?;
        _count(sufr_file, args)
    } else {
        let sufr_file: SufrFile<u64> = SufrFile::read(&args.file)?;
        _count(sufr_file, args)
    }
}

// --------------------------------------------------
fn _count<T>(mut sufr_file: SufrFile<T>, args: &CountArgs) -> Result<()>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + Debug,
{
    let mut output: Box<dyn Write> = match &args.output {
        Some(out_name) => Box::new(File::create(out_name)?),
        _ => Box::new(io::stdout()),
    };

    let queries = parse_locate_queries(&args.query)?;
    let num_queries = queries.len();
    let now = Instant::now();
    let loc_args = SearchOptions {
        queries,
        max_query_len: args.max_query_len,
        low_memory: args.low_memory,
        find_suffixes: false,
    };

    for res in sufr_file.suffix_search(&loc_args)? {
        match res.locations {
            Some(locs) => {
                let ranks = locs.ranks;
                writeln!(output, "{} {}", res.query, ranks.end - ranks.start)?;
            }
            _ => eprintln!("{} not found", res.query),
        }
    }

    info!("Locate of {} finished in {:?}", num_queries, now.elapsed());
    Ok(())
}

// --------------------------------------------------
pub fn create(args: &CreateArgs) -> Result<()> {
    // Read sequence input
    let now = Instant::now();
    let sequence_delimiter = args.sequence_delimiter as u8;
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
fn _create<T>(sa: SufrBuilder<T>, args: &CreateArgs, timer: Instant) -> Result<()>
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
    let text_len = read_text_length(&args.file)? as u64;
    if text_len < u32::MAX as u64 {
        let now = Instant::now();
        let sufr_file: SufrFile<u32> = SufrFile::read(&args.file)?;
        info!("Read sufr file in {:?}", now.elapsed());
        _extract(sufr_file, args)
    } else {
        let now = Instant::now();
        let sufr_file: SufrFile<u64> = SufrFile::read(&args.file)?;
        info!("Read sufr file in {:?}", now.elapsed());
        _extract(sufr_file, args)
    }
}

// --------------------------------------------------
fn _extract<T>(mut sufr_file: SufrFile<T>, args: &ExtractArgs) -> Result<()>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + Debug,
{
    let mut suffixes: Vec<usize> = vec![];
    for suffix in &args.suffixes {
        let mut positions: Vec<_> = parse_pos(suffix)?
            .into_iter()
            .flat_map(|r| r.collect::<Vec<_>>())
            .collect();
        suffixes.append(&mut positions);
    }

    let mut output: Box<dyn Write> = match &args.output {
        Some(out_name) => Box::new(File::create(out_name)?),
        _ => Box::new(io::stdout()),
    };

    let now = Instant::now();
    let text_len = sufr_file.text_len.to_usize();
    let prefix_len = args.prefix_len.unwrap_or(0);
    for suffix in suffixes {
        let start = if suffix - prefix_len > 0 {
            suffix - prefix_len
        } else {
            suffix
        };
        let end = min(
            args.suffix_len.map_or(text_len, |len| suffix + len),
            text_len,
        );
        if let Some(bytes) = sufr_file.text.get(start..end) {
            let lcp = if args.lcp {
                sufr_file
                    .lcp_file
                    .get(suffix)
                    .map_or("NA".to_string(), |v| format!("\t{v}"))
            } else {
                "".to_string()
            };
            let seq = String::from_utf8(bytes.to_vec())?;
            writeln!(output, "{seq}{lcp}")?;
        }
    }
    info!("Extracted suffixes in {:?}", now.elapsed());

    Ok(())
}

// --------------------------------------------------
pub fn list(args: &ListArgs) -> Result<()> {
    let text_len = read_text_length(&args.file)? as u64;
    if text_len < u32::MAX as u64 {
        let now = Instant::now();
        let sa: SufrFile<u32> = SufrFile::read(&args.file)?;
        info!("Read sufr file in {:?}", now.elapsed());
        _list(sa, args)
    } else {
        let now = Instant::now();
        let sa: SufrFile<u64> = SufrFile::read(&args.file)?;
        info!("Read sufr file in {:?}", now.elapsed());
        _list(sa, args)
    }
}

// --------------------------------------------------
fn _list<T>(mut sufr_file: SufrFile<T>, args: &ListArgs) -> Result<()>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + Debug,
{
    let mut ranks: Vec<usize> = vec![];
    for val in &args.ranks {
        let mut parsed: Vec<_> = parse_pos(val)?
            .into_iter()
            .flat_map(|r| r.collect::<Vec<_>>())
            .collect();
        ranks.append(&mut parsed);
    }

    let mut output: Box<dyn Write> = match &args.output {
        Some(out_name) => Box::new(File::create(out_name)?),
        _ => Box::new(io::stdout()),
    };

    let width = sufr_file.text_len.to_string().len();
    let text_len = sufr_file.text_len.to_usize();
    let suffix_len = args.len.unwrap_or(text_len);

    writeln!(output, "{:>width$} {:>width$} {:>width$}", "R", "S", "L")?;

    let mut print = |rank: usize, suffix: usize, lcp: usize| -> Result<()> {
        let end = if suffix + suffix_len > text_len {
            text_len
        } else {
            suffix + suffix_len
        };
        writeln!(
            output,
            "{rank:width$} {suffix:width$} {lcp:width$}: {}",
            String::from_utf8(sufr_file.text[suffix..end].to_vec())?
        )?;
        Ok(())
    };

    if ranks.is_empty() {
        for (rank, (suffix, lcp)) in sufr_file
            .suffix_array_file
            .iter()
            .zip(sufr_file.lcp_file.iter())
            .enumerate()
        {
            print(rank, suffix.to_usize(), lcp.to_usize())?;
        }
    } else {
        for rank in ranks {
            if let (Some(suffix), Some(lcp)) = (
                sufr_file.suffix_array_file.get(rank),
                sufr_file.lcp_file.get(rank),
            ) {
                print(rank, suffix.to_usize(), lcp.to_usize())?;
            } else {
                eprintln!("Invalid rank: {rank}");
            }
        }
    }

    Ok(())
}

// --------------------------------------------------
fn parse_locate_queries(queries: &[String]) -> Result<Vec<String>> {
    let whitespace = Regex::new(r"\s+").unwrap();
    let mut ret = vec![];
    for query in queries {
        if Path::new(&query).exists() {
            let contents = fs::read_to_string(query)?;
            let mut vals: Vec<String> = whitespace
                .split(&contents)
                .filter(|v| !v.is_empty())
                .map(|v| v.to_string())
                .collect();
            ret.append(&mut vals);
        } else {
            ret.push(query.to_string());
        }
    }

    Ok(ret)
}

// --------------------------------------------------
pub fn locate(args: &LocateArgs) -> Result<()> {
    let len = read_text_length(&args.file)? as u64;
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
fn _locate<T>(mut sufr_file: SufrFile<T>, args: &LocateArgs) -> Result<()>
where
    T: Int + FromUsize<T> + Sized + Send + Sync + Debug,
{
    let mut output: Box<dyn Write> = match &args.output {
        Some(out_name) => Box::new(File::create(out_name)?),
        _ => Box::new(io::stdout()),
    };

    let queries = parse_locate_queries(&args.query)?;
    let num_queries = queries.len();
    let now = Instant::now();
    let loc_args = SearchOptions {
        queries,
        max_query_len: args.max_query_len,
        low_memory: args.low_memory,
        find_suffixes: true,
    };

    for mut res in sufr_file.locate(loc_args)? {
        if res.positions.is_empty() {
            eprintln!("{} not found", res.query);
            continue;
        }

        if args.abs {
            writeln!(
                output,
                "{} {}",
                res.query,
                res.positions
                    .into_iter()
                    .map(|p| p.suffix.to_string())
                    .collect::<Vec<_>>()
                    .join(" ")
            )?;
        } else {
            // Sort by name then position
            res.positions.sort_by(|a, b| {
                a.sequence_name
                    .cmp(&b.sequence_name)
                    .then(a.sequence_position.cmp(&b.sequence_position))
            });

            writeln!(output, "{}", res.query)?;
            let mut prev_seq = "".to_string();
            let mut buffer = vec![];
            for pos in res.positions {
                if pos.sequence_name != prev_seq {
                    if !buffer.is_empty() {
                        writeln!(output, "{prev_seq} {}", buffer.join(","))?;
                    }

                    prev_seq = pos.sequence_name;
                    buffer = vec![];
                }
                buffer.push(pos.sequence_position.to_string());
            }

            if !buffer.is_empty() {
                writeln!(output, "{prev_seq} {}", buffer.join(","))?;
            }

            writeln!(output, "//")?;
        }
    }

    info!("Locate of {} finished in {:?}", num_queries, now.elapsed());

    Ok(())
}

// --------------------------------------------------
// Parse an index from a string representation of an integer.
// Ensures the number does not start with '+'.
fn parse_index(input: &str) -> Result<usize> {
    let value_error = || anyhow!(r#"illegal list value: "{input}""#);
    input
        .starts_with('+')
        .then(|| Err(value_error()))
        .unwrap_or_else(|| input.parse::<usize>().map_err(|_| value_error()))
}

// --------------------------------------------------
fn parse_pos(range: &str) -> Result<PositionList> {
    let range_re = Regex::new(r"^(\d+)-(\d+)$").unwrap();
    range
        .split(',')
        .map(|val| {
            parse_index(val).map(|n| n..n + 1).or_else(|e| {
                range_re.captures(val).ok_or(e).and_then(|captures| {
                    let n1 = parse_index(&captures[1])?;
                    let n2 = parse_index(&captures[2])?;
                    if n1 >= n2 {
                        bail!(
                            "First number in range ({n1}) \
                            must be lower than second number ({n2})"
                        )
                    }
                    Ok(n1..n2 + 1)
                })
            })
        })
        .collect::<Result<_, _>>()
        .map_err(From::from)
}

// --------------------------------------------------
pub fn summarize(args: &SummarizeArgs) -> Result<()> {
    let text_len = read_text_length(&args.file)? as u64;
    if text_len < u32::MAX as u64 {
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
fn _summarize<T>(sufr_file: SufrFile<T>, _args: &SummarizeArgs) -> Result<()>
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
        "Allow Ambiguity".to_string(),
        sufr_file.allow_ambiguity.to_string(),
    ]);
    rows.push(vec![
        "Ignore Softmask".to_string(),
        sufr_file.ignore_softmask.to_string(),
    ]);
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

// --------------------------------------------------
#[cfg(test)]
mod tests {
    use super::{parse_index, parse_pos};
    use pretty_assertions::assert_eq;

    #[test]
    fn test_parse_index() {
        let res = parse_index("0");
        assert!(res.is_ok());
        assert_eq!(res.unwrap(), 0);

        let res = parse_index("1");
        assert!(res.is_ok());
        assert_eq!(res.unwrap(), 1);

        let res = parse_index("10");
        assert!(res.is_ok());
        assert_eq!(res.unwrap(), 10);

        let res = parse_index(&usize::MAX.to_string());
        assert!(res.is_ok());
        assert_eq!(res.unwrap(), usize::MAX);

        let res = parse_index("+1");
        assert!(res.is_err());
        assert_eq!(
            res.unwrap_err().to_string(),
            r#"illegal list value: "+1""#.to_string()
        );

        let res = parse_index("-1");
        assert!(res.is_err());
        assert_eq!(
            res.unwrap_err().to_string(),
            r#"illegal list value: "-1""#.to_string()
        );
    }

    #[test]
    fn test_parse_pos() {
        // The empty string is an error
        assert!(parse_pos("").is_err());

        // Zero is OK
        let res = parse_pos("0");
        assert!(res.is_ok());
        assert_eq!(res.unwrap(), vec![0..1]);

        //let res = parse_pos("0-1".to_string());
        //assert!(res.is_err());
        //assert_eq!(res.unwrap_err().to_string(), r#"illegal list value: "0""#);
        //
        //// A leading "+" is an error
        //let res = parse_pos("+1".to_string());
        //assert!(res.is_err());
        //assert_eq!(res.unwrap_err().to_string(), r#"illegal list value: "+1""#,);
        //
        //let res = parse_pos("+1-2".to_string());
        //assert!(res.is_err());
        //assert_eq!(
        //    res.unwrap_err().to_string(),
        //    r#"illegal list value: "+1-2""#,
        //);
        //
        //let res = parse_pos("1-+2".to_string());
        //assert!(res.is_err());
        //assert_eq!(
        //    res.unwrap_err().to_string(),
        //    r#"illegal list value: "1-+2""#,
        //);
        //
        //// Any non-number is an error
        //let res = parse_pos("a".to_string());
        //assert!(res.is_err());
        //assert_eq!(res.unwrap_err().to_string(), r#"illegal list value: "a""#);
        //
        //let res = parse_pos("1,a".to_string());
        //assert!(res.is_err());
        //assert_eq!(res.unwrap_err().to_string(), r#"illegal list value: "a""#);
        //
        //let res = parse_pos("1-a".to_string());
        //assert!(res.is_err());
        //assert_eq!(res.unwrap_err().to_string(), r#"illegal list value: "1-a""#,);
        //
        //let res = parse_pos("a-1".to_string());
        //assert!(res.is_err());
        //assert_eq!(res.unwrap_err().to_string(), r#"illegal list value: "a-1""#,);
        //
        //// Wonky ranges
        //let res = parse_pos("-".to_string());
        //assert!(res.is_err());
        //
        //let res = parse_pos(",".to_string());
        //assert!(res.is_err());
        //
        //let res = parse_pos("1,".to_string());
        //assert!(res.is_err());
        //
        //let res = parse_pos("1-".to_string());
        //assert!(res.is_err());
        //
        //let res = parse_pos("1-1-1".to_string());
        //assert!(res.is_err());
        //
        //let res = parse_pos("1-1-a".to_string());
        //assert!(res.is_err());
        //
        //// First number must be less than second
        //let res = parse_pos("1-1".to_string());
        //assert!(res.is_err());
        //assert_eq!(
        //    res.unwrap_err().to_string(),
        //    "First number in range (1) must be lower than second number (1)"
        //);
        //
        //let res = parse_pos("2-1".to_string());
        //assert!(res.is_err());
        //assert_eq!(
        //    res.unwrap_err().to_string(),
        //    "First number in range (2) must be lower than second number (1)"
        //);
        //
        //// All the following are acceptable
        //let res = parse_pos("1".to_string());
        //assert!(res.is_ok());
        //assert_eq!(res.unwrap(), vec![0..1]);
        //
        //let res = parse_pos("01".to_string());
        //assert!(res.is_ok());
        //assert_eq!(res.unwrap(), vec![0..1]);
        //
        //let res = parse_pos("1,3".to_string());
        //assert!(res.is_ok());
        //assert_eq!(res.unwrap(), vec![0..1, 2..3]);
        //
        //let res = parse_pos("001,0003".to_string());
        //assert!(res.is_ok());
        //assert_eq!(res.unwrap(), vec![0..1, 2..3]);
        //
        //let res = parse_pos("1-3".to_string());
        //assert!(res.is_ok());
        //assert_eq!(res.unwrap(), vec![0..3]);
        //
        //let res = parse_pos("0001-03".to_string());
        //assert!(res.is_ok());
        //assert_eq!(res.unwrap(), vec![0..3]);
        //
        //let res = parse_pos("1,7,3-5".to_string());
        //assert!(res.is_ok());
        //assert_eq!(res.unwrap(), vec![0..1, 6..7, 2..5]);
        //
        //let res = parse_pos("15,19-20".to_string());
        //assert!(res.is_ok());
        //assert_eq!(res.unwrap(), vec![14..15, 18..20]);
    }
}
