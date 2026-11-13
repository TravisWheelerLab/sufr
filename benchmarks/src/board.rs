//! The comparison, drawn as it fills.
//!
//! Every row a run will produce is printed once, empty, and each is rewritten
//! in place as its measurements land. Falls back to one line per finished row
//! where stdout is not a terminal, so piping the run to a file still leaves
//! something to read.

use std::io::{IsTerminal, Write};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::{Duration, Instant};

use crate::process::Timing;

// --- colour

const RESET: &str = "\x1b[0m";
const BOLD: &str = "\x1b[1m";
const DIM: &str = "\x1b[2m";
const GREEN: &str = "\x1b[32m";
const YELLOW: &str = "\x1b[33m";
const RED: &str = "\x1b[31m";
const CYAN: &str = "\x1b[36m";
const MAGENTA: &str = "\x1b[35m";

/// One frame per tick of the painter, so something moves even while a minute
/// long run is in flight.
const SPINNER: [char; 10] = ['⠋', '⠙', '⠹', '⠸', '⠼', '⠴', '⠦', '⠧', '⠇', '⠏'];

/// How often the board repaints itself between events.
const TICK: Duration = Duration::from_millis(80);

const WIDTH: usize = 71;

/// The ratio column's width and its gap, when a group has more than one row.
const RATIO: usize = 9;

/// One colour per binary, so a row is recognisable before it is read.
fn colour_at(n: usize) -> &'static str {
    match n {
        0 => CYAN,
        1 => MAGENTA,
        _ => YELLOW,
    }
}

// --- the rows

/// One row group: an operation and the rows under it, in the order they run.
pub(crate) struct Group {
    pub(crate) op: &'static str,
    pub(crate) why: &'static str,

    /// Each row's label and which binary it belongs to, for its colour.
    pub(crate) rows: Vec<(String, usize)>,
}

struct Row {
    op: &'static str,
    label: String,
    colour: usize,

    /// Every measurement reported so far, one per repeat.
    runs: Vec<Timing>,
}

/// What a finished or partly finished row shows.
struct Stats {
    median_s: f64,
    min_s: f64,
    max_rss_kb: i64,

    /// The median run's, so the two numbers describe one run.
    cpu: f64,
}

impl Row {
    fn stats(&self) -> Option<Stats> {
        if self.runs.is_empty() {
            return None;
        }

        let mut by_wall: Vec<&Timing> = self.runs.iter().collect();
        by_wall.sort_by(|a, b| a.wall_s.total_cmp(&b.wall_s));
        let median = by_wall[by_wall.len() / 2];

        Some(Stats {
            median_s: median.wall_s,
            min_s: by_wall[0].wall_s,
            max_rss_kb: self.runs.iter().map(|t| t.max_rss_kb).max().unwrap_or(0),
            cpu: median.cpu(),
        })
    }
}

/// What is running right now, for the line under the board.
struct Running {
    /// The row on the clock, so it alone spins.
    row: (String, String),
    label: String,
    line: String,
    since: Instant,
}

struct State {
    rows: Vec<Row>,

    /// The command in flight, redrawn every tick so its clock runs.
    running: Option<Running>,

    /// When the whole run started, for the total on the status line.
    started: Instant,

    /// Why each group's rows are doing the same work, printed under the
    /// board because a reader who cannot check the claim cannot read the rows.
    claims: Vec<(&'static str, &'static str)>,
    repeat: usize,
    header: String,

    /// The first binary's label, which the ratio column is against. The
    /// column is there only when a group has more than one row.
    subject: String,
    tty: bool,

    /// Lines the last paint left on screen, which the next one moves back over.
    painted: usize,
}

/// The board, shared with the loop that fills it.
#[derive(Clone)]
pub(crate) struct Board(Arc<Mutex<State>>);

/// The repainting thread, stopped when the run ends.
pub(crate) struct Painter {
    stop: Arc<AtomicBool>,
    handle: Option<thread::JoinHandle<()>>,
}

impl Painter {
    /// Stop repainting and wait for the thread to let go of stdout.
    pub(crate) fn stop(mut self) {
        self.stop.store(true, Ordering::Relaxed);
        if let Some(h) = self.handle.take() {
            let _ = h.join();
        }
    }
}

impl Board {
    pub(crate) fn new(
        plan: &[Group],
        header: String,
        repeat: usize,
        subject: String,
    ) -> Board {
        let mut rows = Vec::new();
        let claims = plan.iter().map(|g| (g.op, g.why)).collect();

        for group in plan {
            for (label, colour) in &group.rows {
                rows.push(Row {
                    op: group.op,
                    label: label.clone(),
                    colour: *colour,
                    runs: Vec::new(),
                });
            }
        }

        Board(Arc::new(Mutex::new(State {
            rows,
            running: None,
            started: Instant::now(),
            claims,
            repeat,
            header,
            subject,
            tty: std::io::stdout().is_terminal(),
            painted: 0,
        })))
    }

    /// Draw the empty board and keep it moving.
    //
    // a thread rather than a repaint per event: the events
    // are seconds apart, and a clock that only ticks when
    // something finishes is not a clock
    pub(crate) fn start(&self) -> Painter {
        let mut state = self.0.lock().unwrap();

        if !state.tty {
            println!("{}", state.header);

            return Painter {
                stop: Arc::new(AtomicBool::new(true)),
                handle: None,
            };
        }

        state.started = Instant::now();
        paint(&mut state);
        drop(state);

        let stop = Arc::new(AtomicBool::new(false));
        let shared = self.0.clone();
        let flag = stop.clone();

        let handle = thread::spawn(move || {
            while !flag.load(Ordering::Relaxed) {
                thread::sleep(TICK);

                if let Ok(mut state) = shared.lock() {
                    paint(&mut state);
                }
            }
        });

        Painter {
            stop,
            handle: Some(handle),
        }
    }

    /// Show `line` as what is on the clock now.
    pub(crate) fn running(&self, op: &str, row: &str, label: String, line: String) {
        let mut state = self.0.lock().unwrap();
        state.running = Some(Running {
            row: (op.to_string(), row.to_string()),
            label,
            line,
            since: Instant::now(),
        });
    }

    /// Add one measurement to a row.
    pub(crate) fn record(&self, op: &str, label: &str, timing: Timing) {
        let mut state = self.0.lock().unwrap();
        let repeat = state.repeat;

        let Some(row) = state
            .rows
            .iter_mut()
            .find(|r| r.op == op && r.label == label)
        else {
            return;
        };

        row.runs.push(timing);
        let full = row.runs.len() == repeat;
        let (op, label) = (row.op, row.label.clone());

        if state.tty {
            paint(&mut state);
        } else if full {
            println!("{}", plain(&line_for(&state, op, &label)));
        }
    }

    /// Leave the finished board in the scrollback, with what each row rests on.
    pub(crate) fn finish(&self) {
        let mut state = self.0.lock().unwrap();
        state.running = None;

        if state.tty {
            paint(&mut state);
        }

        println!();
        for (op, why) in &state.claims {
            let line = format!("  {DIM}{op:<14} {why}{RESET}");
            println!("{}", if state.tty { line } else { plain(&line) });
        }
        println!();
    }
}

// --- drawing

fn line_for(state: &State, op: &str, label: &str) -> String {
    let row = state
        .rows
        .iter()
        .find(|r| r.op == op && r.label == label)
        .expect("a row that was just written");

    let first = state
        .rows
        .iter()
        .find(|r| r.op == op)
        .is_some_and(|r| std::ptr::eq(r, row));

    let group = if first { op } else { "" };
    let colour = colour_at(row.colour);
    let weight = if first { BOLD } else { "" };

    let Some(stats) = row.stats() else {
        return format!("  {DIM}{group:<14} {label:<17}          —{RESET}");
    };

    if row.runs.len() < state.repeat {
        let on_clock = state
            .running
            .as_ref()
            .is_some_and(|r| r.row.0 == op && r.row.1 == label);

        // only the row on the clock spins; a row that is
        // partly done and waiting shows its count
        let frame = if on_clock {
            SPINNER[(state.started.elapsed().as_secs_f64() / TICK.as_secs_f64())
                as usize
                % SPINNER.len()]
        } else {
            ' '
        };

        return format!(
            "  {group:<14} {weight}{colour}{label:<17}{RESET} {colour}{frame}{RESET} \
             {DIM}{}/{}{RESET}",
            row.runs.len(),
            state.repeat,
        );
    }

    // the first row of the group is what every ratio in it
    // is against: the question a group answers is what the
    // other binaries cost relative to the one written first
    let ratio = if compared(state) {
        let subject = state
            .rows
            .iter()
            .find(|r| r.op == op)
            .and_then(|r| r.stats())
            .map(|s| stats.median_s / s.median_s);

        match subject {
            Some(r) if (r - 1.0).abs() < 0.005 => {
                format!(" {DIM}{:>8}{RESET}", "1.00x")
            }
            Some(r) if r > 1.0 => format!(" {RED}{:>8}{RESET}", format!("{r:.2}x")),
            Some(r) => format!(" {GREEN}{:>8}{RESET}", format!("{r:.2}x")),
            None => format!(" {:>8}", "—"),
        }
    } else {
        String::new()
    };

    format!(
        "  {group:<14} {weight}{colour}{label:<17} {:>10}{RESET} {DIM}{:>10}{RESET}{weight}{colour} {:>10} {:>5}{RESET}{ratio}",
        seconds(stats.median_s),
        seconds(stats.min_s),
        // ru_maxrss is KiB on linux, so this is MiB and
        // says so: the difference reaches 5% by a gigabyte
        format!("{:.0} MiB", stats.max_rss_kb as f64 / 1024.0),
        format!("{:.0}%", stats.cpu * 100.0),
    )
}

/// Whether any group has more than one row.
fn compared(state: &State) -> bool {
    state
        .rows
        .iter()
        .any(|r| state.rows.iter().filter(|o| o.op == r.op).count() > 1)
}

fn seconds(s: f64) -> String {
    format!("{s:.2} s")
}

/// Rewrite the whole board in place.
//
// the cursor goes back over what was painted and every line
// is cleared before it is written, because a shorter line
// would otherwise leave the tail of a longer one behind
fn paint(state: &mut State) {
    let mut out = std::io::stdout().lock();

    if state.painted > 0 {
        let _ = write!(out, "\x1b[{}A", state.painted);
    }

    let (ratio, width) = if compared(state) {
        (
            format!(" {:>8}", format!("vs {}", state.subject)),
            WIDTH + RATIO,
        )
    } else {
        (String::new(), WIDTH)
    };

    let mut lines = vec![
        format!("\x1b[2K {DIM}{}{RESET}", state.header),
        "\x1b[2K".to_string(),
        format!(
            "\x1b[2K  {BOLD}{:<14} {:<17} {:>10} {:>10} {:>10} {:>5}{ratio}{RESET}",
            "operation", "row", "median", "min", "rss", "cpu"
        ),
        format!("\x1b[2K  {DIM}{}{RESET}", "─".repeat(width)),
    ];

    let ops: Vec<&str> = {
        let mut seen = Vec::new();
        for r in &state.rows {
            if !seen.contains(&r.op) {
                seen.push(r.op);
            }
        }
        seen
    };

    for (n, op) in ops.iter().enumerate() {
        // a rule between groups, so the eye reads one
        // operation's rows together rather than a column
        if n > 0 {
            lines.push(format!("\x1b[2K  {DIM}{}{RESET}", "┄".repeat(width)));
        }

        for label in state
            .rows
            .iter()
            .filter(|r| r.op == *op)
            .map(|r| r.label.clone())
            .collect::<Vec<_>>()
        {
            lines.push(format!("\x1b[2K{}", line_for(state, op, &label)));
        }
    }

    lines.push("\x1b[2K".to_string());
    for line in status(state) {
        lines.push(format!("\x1b[2K{line}"));
    }

    for line in &lines {
        let _ = writeln!(out, "{line}");
    }
    let _ = out.flush();

    state.painted = lines.len();
}

/// The line under the board: what is running, how long it has been, and how
/// far through the whole run this is.
fn status(state: &State) -> Vec<String> {
    let done = state.rows.iter().map(|r| r.runs.len()).sum::<usize>();
    let want = state.rows.len() * state.repeat;
    let total = state.started.elapsed().as_secs_f64();

    let counter = format!("  {DIM}{done}/{want} runs · {total:.1}s{RESET}");

    let Some(running) = &state.running else {
        return vec![counter];
    };

    let frame = SPINNER[(total / TICK.as_secs_f64()) as usize % SPINNER.len()];
    let held = running.since.elapsed().as_secs_f64();

    let shown = short(&running.line);

    vec![
        format!(
            "  {MAGENTA}{frame}{RESET} {BOLD}{:<24}{RESET}{DIM}{held:>6.1}s{RESET}  {CYAN}{shown}{RESET}",
            running.label,
        ),
        counter,
    ]
}

/// A command as something worth watching go past: every path down to its last
/// segment.
fn short(line: &str) -> String {
    let words: Vec<String> = line
        .split_whitespace()
        .map(|word| match word.rsplit_once('/') {
            Some((_, last)) if !last.is_empty() => last.to_string(),
            _ => word.to_string(),
        })
        .collect();

    let shown = words.join(" ");
    let width = 76;

    if shown.chars().count() > width {
        format!("{}…", shown.chars().take(width).collect::<String>())
    } else {
        shown
    }
}

/// A line with its colour taken back out, for a run that is not going to a
/// terminal.
//
// stripped after formatting rather than never added:
// padding is applied to the values, so the columns line up
// either way and there is one code path to read
fn plain(line: &str) -> String {
    let mut out = String::with_capacity(line.len());
    let mut chars = line.chars();

    while let Some(c) = chars.next() {
        if c != '\x1b' {
            out.push(c);
            continue;
        }

        // an escape runs to the letter that ends it
        for c in chars.by_ref() {
            if c.is_ascii_alphabetic() {
                break;
            }
        }
    }

    out
}
