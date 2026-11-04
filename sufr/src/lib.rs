use anyhow::{anyhow, bail, Result};
use clap::{builder::PossibleValue, Parser, ValueEnum};
use format_num::NumberFormat;
use libsufr::{
    read_queries, Alphabet, Context, Hits, OpenOptions, SeedMask, Storage, SuffixArray,
    SuffixArrayBuilder, SuffixSortType,
};
use log::info;
use rayon::prelude::*;
use regex::Regex;
use std::{
    ffi::OsStr,
    fmt::Debug,
    fs::{self, File},
    io::{self, Write},
    ops::Range,
    path::PathBuf,
    str,
    time::Instant,
};
use tabled::Table;

type PositionList = Vec<Range<usize>>;

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
#[command(about, alias = "cr")]
pub struct CreateArgs {
    /// Input file
    #[arg(value_name = "INPUT")]
    pub input: String,

    /// Subproblem count
    #[arg(short, long, value_name = "NUM_PARTS", default_value = "1024")]
    pub num_partitions: usize,

    /// Max context
    #[arg(short, long, value_name = "CONTEXT", conflicts_with = "seed_mask")]
    pub max_query_len: Option<usize>,

    /// Output file
    #[arg(short, long, value_name = "OUTPUT")]
    pub output: Option<String>,

    /// Input is DNA
    #[arg(short('d'), long("dna"))]
    pub is_dna: bool,

    /// Allow suffixes starting with ambiguity codes
    #[arg(short, long, requires = "is_dna")]
    pub allow_ambiguity: bool,

    /// Ignore suffixes in soft-mask/lowercase regions
    #[arg(short, long, requires = "is_dna")]
    pub ignore_softmask: bool,

    /// Character to separate sequences
    #[arg(short('D'), long, default_value = "%", value_name = "DELIM")]
    pub sequence_delimiter: char,

    /// Spaced seeds mask
    #[arg(short, long, value_name = "MASK")]
    pub seed_mask: Option<String>,

    /// Also write the LCP array (needed by "list --show-lcp"; speeds up
    /// query-time max query lengths shorter than the build-time value)
    #[arg(long)]
    pub write_lcp: bool,
}

#[derive(Debug, Parser)]
#[command(about, alias = "ex")]
pub struct ExtractArgs {
    /// Maximum query length
    #[arg(short, long, value_name = "LEN")]
    pub max_query_len: Option<usize>,

    /// Low memory
    #[arg(short, long)]
    pub low_memory: bool,

    /// Very low memory
    #[arg(short, long, conflicts_with = "low_memory")]
    pub very_low_memory: bool,

    /// Prefix length
    #[arg(short, long, value_name = "PREFIX_LEN")]
    pub prefix_len: Option<usize>,

    /// Suffix length
    #[arg(short, long, value_name = "SUFFIX_LEN")]
    pub suffix_len: Option<usize>,

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
#[command(about, alias = "ls")]
pub struct ListArgs {
    /// Sufr file
    #[arg(value_name = "FILE")]
    pub file: String,

    /// Ranks of suffixes to show
    #[arg(value_name = "RANK")]
    pub ranks: Vec<String>,

    /// Show rank column
    #[arg(short('r'), long)]
    pub show_rank: bool,

    /// Show suffix position column
    #[arg(short('s'), long)]
    pub show_suffix: bool,

    /// Show LCP column
    #[arg(short('p'), long)]
    pub show_lcp: bool,

    /// Very low memory
    #[arg(short, long)]
    pub very_low_memory: bool,

    /// Length of suffixes to show
    #[arg(long, value_name = "LEN")]
    pub len: Option<usize>,

    /// Number of suffixes to show
    #[arg(short, long, value_name = "LEN")]
    pub number: Option<usize>,

    /// Output
    #[arg(short, long, value_name = "OUT")]
    pub output: Option<String>,
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

    /// Low memory
    #[arg(short, long)]
    pub low_memory: bool,

    /// Very low memory
    #[arg(short, long, conflicts_with = "low_memory")]
    pub very_low_memory: bool,

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
    /// Output
    #[arg(short, long, value_name = "OUT")]
    pub output: Option<String>,

    /// Maximum query length
    #[arg(short, long, value_name = "LEN")]
    pub max_query_len: Option<usize>,

    /// Low memory
    #[arg(short, long)]
    pub low_memory: bool,

    /// Very low memory
    #[arg(short, long, conflicts_with = "low_memory")]
    pub very_low_memory: bool,

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

/// Open a _.sufr_ file for `count`, `extract` or `locate`.
fn open_for_search(
    file: &str,
    low_memory: bool,
    very_low_memory: bool,
    max_query_len: Option<usize>,
) -> Result<SuffixArray> {
    let text = if very_low_memory {
        Storage::Disk
    } else {
        Storage::Memory
    };
    let suffixes = if low_memory || very_low_memory {
        Storage::Disk
    } else {
        Storage::Memory
    };
    SuffixArray::open_with(
        file,
        &OpenOptions {
            text,
            suffixes,
            max_query_len,
        },
    )
}

/// The named output file, or stdout.
fn create_output(output: &Option<String>) -> Result<Box<dyn Write>> {
    match output {
        Some(out_name) => Ok(Box::new(File::create(out_name)?)),
        _ => Ok(Box::new(io::stdout())),
    }
}

pub fn count(args: &CountArgs) -> Result<()> {
    let suffix_array = open_for_search(
        &args.file,
        args.low_memory,
        args.very_low_memory,
        args.max_query_len,
    )?;
    let mut output = create_output(&args.output)?;
    let queries = read_queries(&args.query)?;
    let now = Instant::now();

    let hits = suffix_array.par_search_all(&queries)?;
    for (query, hits) in queries.iter().zip(hits) {
        writeln!(output, "{query} {}", hits.count())?;
    }

    info!("Count of {} finished in {:?}", queries.len(), now.elapsed());
    Ok(())
}

pub fn create(args: &CreateArgs) -> Result<()> {
    let outfile = &args.output.clone().unwrap_or(format!(
        "{}.sufr",
        PathBuf::from(&args.input)
            .file_stem()
            .unwrap_or(OsStr::new("out"))
            .to_string_lossy()
    ));

    let alphabet = if args.is_dna {
        Alphabet::Dna {
            allow_ambiguity: args.allow_ambiguity,
            ignore_softmask: args.ignore_softmask,
        }
    } else {
        Alphabet::Bytes
    };
    let mut builder = SuffixArrayBuilder::from_sequence_file(
        &args.input,
        args.sequence_delimiter as u8,
    )?
    .with_alphabet(alphabet)
    .with_num_partitions(args.num_partitions);

    match (args.max_query_len, &args.seed_mask) {
        (Some(_), Some(_)) => bail!("Cannot use max_query_len and seed_mask together"),
        (Some(n), None) => {
            builder = builder.with_sort_type(SuffixSortType::MaxQueryLen(n))
        }
        (None, Some(mask)) => {
            builder = builder.with_sort_type(SuffixSortType::Mask(SeedMask::new(mask)?))
        }
        (None, None) => (),
    }

    builder = builder.with_write_lcp(args.write_lcp);

    let now = Instant::now();
    builder.build_to_file(outfile, true)?;

    let meta = fs::metadata(outfile)?;
    let bytes_written = meta.len();

    let num_fmt = NumberFormat::new();
    info!(
        "Wrote {} byte{} to '{outfile}' in {:?}",
        num_fmt.format(",.0", bytes_written as f64),
        if bytes_written == 1 { "" } else { "s" },
        now.elapsed()
    );

    Ok(())
}

pub fn extract(args: &ExtractArgs) -> Result<()> {
    let suffix_array = open_for_search(
        &args.file,
        args.low_memory,
        args.very_low_memory,
        args.max_query_len,
    )?;
    let queries = read_queries(&args.query)?;
    let now = Instant::now();
    let mut output = create_output(&args.output)?;

    let hits = suffix_array.par_search_all(&queries)?;
    for (query, hits) in queries.iter().zip(hits) {
        let context = Context {
            before: args.prefix_len.unwrap_or(0),
            after: args.suffix_len,
        };
        let extractions = hits.extract(context)?;
        if extractions.is_empty() {
            eprintln!("{query} not found");
        }
        for ex in extractions {
            let start = ex.location.sequence_start;
            let window = ex.sequence_range;
            let text = suffix_array.text(start + window.start..start + window.end)?;
            writeln!(
                output,
                ">{}:{}-{} {query} {}\n{}",
                ex.location.sequence_name,
                window.start,
                window.end,
                ex.suffix_offset,
                str::from_utf8(&text)?
            )?;
        }
    }
    info!("Extract finished in {:?}", now.elapsed());

    Ok(())
}

pub fn list(args: &ListArgs) -> Result<()> {
    let text = if args.very_low_memory {
        Storage::Disk
    } else {
        Storage::Memory
    };
    let suffix_array = SuffixArray::open_with(
        &args.file,
        &OpenOptions {
            text,
            ..OpenOptions::default()
        },
    )?;
    let ranks = parse_ranks(&args.ranks)?;

    if args.show_lcp && suffix_array.header().lcp_pos.is_none() {
        bail!(
            "{} has no LCP array; rebuild with --write-lcp to show LCP values",
            args.file
        );
    }

    // TODO: list -o is accepted and ignored
    let mut output = io::stdout();
    let text_len = suffix_array.header().text_len;
    let num_suffixes = suffix_array.header().num_suffixes;
    let width = text_len.to_string().len();
    let row = |rank: usize, suffix: usize, lcp: Option<usize>| -> Result<String> {
        let mut row = String::new();
        if args.show_rank {
            row += &format!("{rank:width$} ");
        }
        if args.show_suffix {
            row += &format!("{suffix:width$} ");
        }
        if let Some(lcp) = lcp {
            row += &format!("{lcp:width$} ");
        }
        let end = args.len.map_or(text_len, |len| suffix + len);
        row += &String::from_utf8(suffix_array.text(suffix..end)?)?;
        Ok(row)
    };

    if ranks.is_empty() {
        let end = match args.number {
            Some(number) if number > 0 => number.min(num_suffixes),
            _ => num_suffixes,
        };

        // a chunk of ranks at a time, so a file of billions of
        // suffixes is never read whole
        const CHUNK: usize = 1 << 16;
        for start in (0..end).step_by(CHUNK) {
            let ranks = start..(start + CHUNK).min(end);
            let suffixes = suffix_array.suffixes(ranks.clone())?;
            let lcps = if args.show_lcp {
                suffix_array.lcps(ranks.clone())?
            } else {
                vec![]
            };
            for (i, rank) in ranks.enumerate() {
                writeln!(output, "{}", row(rank, suffixes[i], lcps.get(i).copied())?)?;
            }
        }
    } else {
        for rank in ranks {
            if rank >= num_suffixes {
                eprintln!("Invalid rank: {rank}");
                continue;
            }
            let suffix = suffix_array.suffixes(rank..rank + 1)?[0];
            let lcp = if args.show_lcp {
                Some(suffix_array.lcps(rank..rank + 1)?[0])
            } else {
                None
            };
            writeln!(output, "{}", row(rank, suffix, lcp)?)?;
        }
    }

    Ok(())
}

pub fn locate(args: &LocateArgs) -> Result<()> {
    let suffix_array = open_for_search(
        &args.file,
        args.low_memory,
        args.very_low_memory,
        args.max_query_len,
    )?;
    let mut output = create_output(&args.output)?;
    let queries = read_queries(&args.query)?;
    let now = Instant::now();

    let hits = suffix_array.par_search_all(&queries)?;
    if args.abs {
        write_suffixes(&mut output, &queries, &hits)?;
    } else {
        write_locations(&mut output, &queries, &hits)?;
    }

    info!(
        "Locate of {} finished in {:?}",
        queries.len(),
        now.elapsed()
    );

    Ok(())
}

/// Write each query's hits as positions in the text.
fn write_suffixes(
    output: &mut dyn Write,
    queries: &[String],
    hits: &[Hits],
) -> Result<()> {
    let suffixes: Vec<_> = hits.par_iter().map(|hits| hits.suffixes()).collect();
    for (query, suffixes) in queries.iter().zip(suffixes) {
        let suffixes = suffixes?;
        if suffixes.is_empty() {
            eprintln!("{query} not found");
            continue;
        }
        let suffixes: Vec<_> = suffixes.iter().map(|s| s.to_string()).collect();
        writeln!(output, "{query} {}", suffixes.join(" "))?;
    }
    Ok(())
}

/// Write each query's hits as positions within their sequences.
fn write_locations(
    output: &mut dyn Write,
    queries: &[String],
    hits: &[Hits],
) -> Result<()> {
    let locations: Vec<_> = hits.par_iter().map(|hits| hits.locate()).collect();
    for (query, locations) in queries.iter().zip(locations) {
        let locations = locations?;
        if locations.is_empty() {
            eprintln!("{query} not found");
            continue;
        }
        writeln!(output, "{query}")?;
        for sequence in locations.chunk_by(|a, b| a.sequence_name == b.sequence_name) {
            let positions: Vec<_> = sequence
                .iter()
                .map(|loc| loc.position.to_string())
                .collect();
            writeln!(
                output,
                "{} {}",
                sequence[0].sequence_name,
                positions.join(",")
            )?;
        }
        writeln!(output, "//")?;
    }
    Ok(())
}

/// The ranks named by each of `ranks`, in order.
fn parse_ranks(ranks: &[String]) -> Result<Vec<usize>> {
    let mut parsed = vec![];
    for val in ranks {
        for range in parse_pos(val)? {
            parsed.extend(range);
        }
    }
    Ok(parsed)
}

/// Parse a rank that does not start with `+`.
fn parse_index(input: &str) -> Result<usize> {
    let value_error = || anyhow!(r#"illegal list value: "{input}""#);
    if input.starts_with('+') {
        Err(value_error())
    } else {
        input.parse::<usize>().map_err(|_| value_error())
    }
}

/// Parse a comma-separated list of ranks and `first-last` rank ranges.
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
}

pub fn summarize(args: &SummarizeArgs) -> Result<()> {
    let suffix_array = SuffixArray::open(&args.file)?;
    let rows = suffix_array.metadata()?.summary();
    let table = Table::from_iter(
        rows.into_iter()
            .map(|(label, value)| [label.to_string(), value]),
    );
    println!("{table}");

    Ok(())
}

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
    }
}
