//! One measured process, and the page cache around it.

use std::fs::File;
use std::hash::{DefaultHasher, Hasher};
use std::io::Read;
use std::os::fd::AsRawFd;
use std::path::Path;
use std::process::Command;
use std::time::Instant;

use anyhow::{Context, bail};

/// What one run of a command cost.
#[derive(Clone, Copy)]
pub(crate) struct Timing {
    pub(crate) wall_s: f64,
    pub(crate) user_s: f64,
    pub(crate) sys_s: f64,
    pub(crate) max_rss_kb: i64,
}

impl Timing {
    /// Cores kept busy over the run: (user + sys) / wall.
    pub(crate) fn cpu(&self) -> f64 {
        (self.user_s + self.sys_s) / self.wall_s
    }
}

/// Run `cmd` to completion and report what it cost.
pub(crate) fn measure(cmd: &mut Command) -> anyhow::Result<Timing> {
    let started = Instant::now();
    let child = cmd
        .spawn()
        .with_context(|| format!("spawning {}", cmd.get_program().to_string_lossy()))?;
    let pid = child.id() as libc::pid_t;

    let mut status: libc::c_int = 0;
    // SAFETY: an all-zero rusage is a valid value, and both
    //         locals outlive the call that writes them
    let mut usage: libc::rusage = unsafe { std::mem::zeroed() };
    let waited = unsafe { libc::wait4(pid, &mut status, 0, &mut usage) };
    let wall_s = started.elapsed().as_secs_f64();

    if waited != pid {
        bail!("wait4: {}", std::io::Error::last_os_error());
    }
    if !libc::WIFEXITED(status) || libc::WEXITSTATUS(status) != 0 {
        bail!(
            "{} exited with {status}; its stderr is in runs/stderr",
            cmd.get_program().to_string_lossy()
        );
    }

    Ok(Timing {
        wall_s,
        user_s: seconds(usage.ru_utime),
        sys_s: seconds(usage.ru_stime),
        max_rss_kb: usage.ru_maxrss,
    })
}

fn seconds(t: libc::timeval) -> f64 {
    t.tv_sec as f64 + t.tv_usec as f64 / 1e6
}

/// Drop each file's pages from the page cache, so the next read of it reaches
/// the disk.
pub(crate) fn evict(paths: &[impl AsRef<Path>]) -> anyhow::Result<()> {
    for path in paths {
        let path = path.as_ref();
        let Ok(file) = File::open(path) else {
            continue;
        };

        // POSIX_FADV_DONTNEED drops clean pages and leaves
        // dirty ones, so a file written moments ago would
        // survive its own eviction without this
        file.sync_all()?;

        // SAFETY: the fd is open for the length of the call,
        //         and posix_fadvise reads no memory through
        //         it; 0, 0 names the whole file
        let rc = unsafe {
            libc::posix_fadvise(file.as_raw_fd(), 0, 0, libc::POSIX_FADV_DONTNEED)
        };
        if rc != 0 {
            bail!("posix_fadvise({}) failed with {rc}", path.display());
        }
    }

    Ok(())
}

/// One fingerprint over the bytes of every file in `paths` that exists.
pub(crate) fn fingerprint(paths: &[impl AsRef<Path>]) -> anyhow::Result<u64> {
    let mut hasher = DefaultHasher::new();
    let mut buffer = vec![0u8; 1 << 20];

    for path in paths {
        let path = path.as_ref();
        let Ok(mut file) = File::open(path) else {
            continue;
        };

        loop {
            let n = file
                .read(&mut buffer)
                .with_context(|| format!("reading {}", path.display()))?;
            if n == 0 {
                break;
            }
            hasher.write(&buffer[..n]);
        }
    }

    Ok(hasher.finish())
}
