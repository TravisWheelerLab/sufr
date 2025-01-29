use anyhow::{anyhow, bail, Result};
use clap::{builder::PossibleValue, Parser, ValueEnum};
use format_num::NumberFormat;
use libsufr::{
    suffix_array::SuffixArray,
    types::{
        CountOptions, ExtractOptions, ListOptions, LocateOptions, SuffixSortType,
        SufrBuilderArgs,
    },
    util::read_sequence_file,
};
use log::info;
use regex::Regex;
use std::{
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

    // Check sufr file for correctness
    //Check(CheckArgs),
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

//#[derive(Debug, Parser)]
//#[command(about, alias = "ch")]
//pub struct CheckArgs {
//    /// Sufr file
//    #[arg(value_name = "SUFR")]
//    pub file: String,
//
//    /// List errors
//    #[arg(short, long)]
//    pub verbose: bool,
//}

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
    #[arg(short, long, value_name = "CONTEXT", conflicts_with = "seed_mask")]
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

    /// Spaced seeds mask
    #[arg(short, long, value_name = "MASK")]
    pub seed_mask: Option<String>,

    /// Random seed
    #[arg(short, long, value_name = "RANDSEED", default_value = "42")]
    pub random_seed: u64,
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
pub fn count(args: &CountArgs) -> Result<()> {
    let mut suffix_array = SuffixArray::read(&args.file, args.very_low_memory)?;
    let mut output: Box<dyn Write> = match &args.output {
        Some(out_name) => Box::new(File::create(out_name)?),
        _ => Box::new(io::stdout()),
    };

    let queries = parse_locate_queries(&args.query)?;
    let num_queries = queries.len();
    let now = Instant::now();
    let count_args = CountOptions {
        queries,
        max_query_len: args.max_query_len,
        low_memory: if args.very_low_memory {
            true
        } else {
            args.low_memory
        },
    };

    for res in suffix_array.count(count_args)? {
        writeln!(output, "{} {}", res.query, res.count)?;
    }

    info!("Count of {} finished in {:?}", num_queries, now.elapsed());
    Ok(())
}

// --------------------------------------------------
pub fn create(args: &CreateArgs) -> Result<()> {
    // Read sequence input
    let now = Instant::now();
    let sequence_delimiter = args.sequence_delimiter as u8;
    let seq_data = read_sequence_file(Path::new(&args.input), sequence_delimiter)?;
    let text_len = seq_data.seq.len();
    let num_fmt = NumberFormat::new();
    info!(
        "Read input of len {} in {:?}",
        num_fmt.format(",.0", text_len as f64),
        now.elapsed()
    );

    let outfile = &args.output.clone().unwrap_or(format!(
        "{}.sufr",
        PathBuf::from(&args.input)
            .file_stem()
            .unwrap_or(OsStr::new("out"))
            .to_string_lossy()
    ));

    let builder_args = SufrBuilderArgs {
        text: seq_data.seq,
        path: Some(outfile.to_string()),
        low_memory: true,
        max_query_len: args.max_query_len,
        is_dna: args.is_dna,
        allow_ambiguity: args.allow_ambiguity,
        ignore_softmask: args.ignore_softmask,
        sequence_starts: seq_data.start_positions.into_iter().collect(),
        sequence_names: seq_data.sequence_names,
        num_partitions: args.num_partitions,
        seed_mask: args.seed_mask.clone(),
        random_seed: args.random_seed,
    };

    let now = Instant::now();
    let path = SuffixArray::write(builder_args)?;
    let meta = fs::metadata(&path)?;
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

// --------------------------------------------------
pub fn extract(args: &ExtractArgs) -> Result<()> {
    let mut suffix_array = SuffixArray::read(&args.file, args.very_low_memory)?;
    let queries = parse_locate_queries(&args.query)?;
    let now = Instant::now();
    let extract_args = ExtractOptions {
        queries,
        max_query_len: args.max_query_len,
        low_memory: if args.very_low_memory {
            true
        } else {
            args.low_memory
        },
        prefix_len: args.prefix_len,
        suffix_len: args.suffix_len,
    };

    let mut output: Box<dyn Write> = match &args.output {
        Some(out_name) => Box::new(File::create(out_name)?),
        _ => Box::new(io::stdout()),
    };
    for res in suffix_array.extract(extract_args)? {
        if res.sequences.is_empty() {
            eprintln!("{} not found", res.query);
        } else {
            for seq in res.sequences {
                writeln!(
                    output,
                    ">{}:{}-{} {} {}\n{}",
                    seq.sequence_name,
                    seq.sequence_range.start,
                    seq.sequence_range.end,
                    res.query,
                    seq.suffix_offset,
                    suffix_array.string_at(
                        seq.sequence_start + seq.sequence_range.start,
                        Some(seq.sequence_range.end - seq.sequence_range.start),
                    )?
                )?;
            }
        }
    }
    info!("Extract finished in {:?}", now.elapsed());

    Ok(())
}

// --------------------------------------------------
pub fn list(args: &ListArgs) -> Result<()> {
    let mut suffix_array = SuffixArray::read(&args.file, args.very_low_memory)?;
    let mut ranks: Vec<usize> = vec![];
    for val in &args.ranks {
        let mut parsed: Vec<_> = parse_pos(val)?
            .into_iter()
            .flat_map(|r| r.collect::<Vec<_>>())
            .collect();
        ranks.append(&mut parsed);
    }

    let list_opts = ListOptions {
        ranks,
        show_rank: args.show_rank,
        show_suffix: args.show_suffix,
        show_lcp: args.show_lcp,
        len: args.len,
        number: args.number,
        // defaults to stdout
        output: None,
    };
    suffix_array.list(list_opts)?;

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
    let mut suffix_array = SuffixArray::read(&args.file, args.very_low_memory)?;
    let mut output: Box<dyn Write> = match &args.output {
        Some(out_name) => Box::new(File::create(out_name)?),
        _ => Box::new(io::stdout()),
    };

    let queries = parse_locate_queries(&args.query)?;
    let num_queries = queries.len();
    let now = Instant::now();
    let loc_args = LocateOptions {
        queries,
        max_query_len: args.max_query_len,
        low_memory: if args.very_low_memory {
            true
        } else {
            args.low_memory
        },
    };

    for mut res in suffix_array.locate(loc_args)? {
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
    let suffix_array = SuffixArray::read(&args.file, true)?;
    let meta = suffix_array.metadata()?;
    let num_fmt = NumberFormat::new();
    let mut rows = vec![vec!["Filename".to_string(), meta.filename.clone()]];
    rows.push(vec![
        "Modified".to_string(),
        meta.modified.format("%Y-%m-%d %H:%M").to_string(),
    ]);
    rows.push(vec![
        "File Size".to_string(),
        format!("{} bytes", num_fmt.format(",.0", meta.file_size as f64)),
    ]);
    rows.push(vec![
        "File Version".to_string(),
        meta.file_version.to_string(),
    ]);
    rows.push(vec!["DNA".to_string(), meta.is_dna.to_string()]);
    rows.push(vec![
        "Allow Ambiguity".to_string(),
        meta.allow_ambiguity.to_string(),
    ]);
    rows.push(vec![
        "Ignore Softmask".to_string(),
        meta.ignore_softmask.to_string(),
    ]);
    rows.push(vec![
        "Text Length".to_string(),
        num_fmt.format(",.0", meta.text_len as f64),
    ]);
    rows.push(vec![
        "Len Suffixes".to_string(),
        num_fmt.format(",.0", meta.len_suffixes as f64),
    ]);

    match meta.sort_type {
        SuffixSortType::Mask(seed_mask) => {
            rows.push(vec!["Seed mask".to_string(), seed_mask.mask])
        }
        SuffixSortType::MaxQueryLen(max_query_len) => rows.push(vec![
            "Max query len".to_string(),
            num_fmt.format(",.0", max_query_len as f64),
        ]),
    };

    rows.push(vec![
        "Num sequences".to_string(),
        num_fmt.format(",.0", meta.num_sequences as f64),
    ]);
    let seq_starts = meta
        .sequence_starts
        .into_iter()
        .map(|v| v.to_string())
        .collect::<Vec<_>>()
        .join(", ");
    rows.push(vec![
        "Sequence starts".to_string(),
        textwrap::wrap(&seq_starts, 40).join("\n"),
    ]);
    let seq_names = meta.sequence_names.join(", ");
    rows.push(vec![
        "Sequence names".to_string(),
        textwrap::wrap(&seq_names, 40).join("\n"),
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
