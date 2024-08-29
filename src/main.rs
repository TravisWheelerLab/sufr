use anyhow::Result;
use clap::Parser;
use sufr::{Cli, Command};

// --------------------------------------------------
fn main() {
    if let Err(e) = run(Cli::parse()) {
        eprintln!("Error: {e}");
        std::process::exit(1);
    }
}

// --------------------------------------------------
fn run(args: Cli) -> Result<()> {
    match &args.command {
        Some(Command::Check(args)) => {
            sufr::check(args)?;
            Ok(())
        }
        Some(Command::Create(args)) => {
            sufr::create(args)?;
            Ok(())
        }
        Some(Command::Read(args)) => {
            sufr::read(args)?;
            Ok(())
        }
        _ => unreachable!(),
    }
}
