//! The memory available to this process

/// The memory this process can allocate without swapping, in bytes, or
/// `None` where it is not known (only Linux is supported).
pub(crate) fn available_memory() -> Option<u64> {
    #[cfg(target_os = "linux")]
    {
        linux::available_memory()
    }
    #[cfg(not(target_os = "linux"))]
    {
        None
    }
}

#[cfg(target_os = "linux")]
mod linux {
    use std::{
        fs,
        path::{Path, PathBuf},
    };

    pub(super) fn available_memory() -> Option<u64> {
        match (mem_available(), cgroup_headroom()) {
            (Some(system), Some(cgroup)) => Some(system.min(cgroup)),
            (system, cgroup) => system.or(cgroup),
        }
    }

    /// The system's `MemAvailable`, in bytes.
    fn mem_available() -> Option<u64> {
        let meminfo = fs::read_to_string("/proc/meminfo").ok()?;
        let kib: u64 = meminfo
            .lines()
            .find_map(|line| line.strip_prefix("MemAvailable:"))?
            .split_whitespace()
            .next()?
            .parse()
            .ok()?;
        Some(kib * 1024)
    }

    /// The smallest gap between a memory limit and its usage among this
    /// process's cgroups and their ancestors, in bytes.
    //
    // batch schedulers such as slurm set the job's limit on
    // an ancestor of the cgroup the process runs in, and
    // MemAvailable does not reflect it
    fn cgroup_headroom() -> Option<u64> {
        let cgroups = fs::read_to_string("/proc/self/cgroup").ok()?;
        let mut headroom: Option<u64> = None;
        for line in cgroups.lines() {
            // "hierarchy-id:controllers:path", where the
            // controllers are empty for cgroup v2
            let mut fields = line.splitn(3, ':');
            let (Some(_), Some(controllers), Some(path)) =
                (fields.next(), fields.next(), fields.next())
            else {
                continue;
            };
            let (root, limit_file, usage_file) = if controllers.is_empty() {
                ("/sys/fs/cgroup", "memory.max", "memory.current")
            } else if controllers.split(',').any(|c| c == "memory") {
                (
                    "/sys/fs/cgroup/memory",
                    "memory.limit_in_bytes",
                    "memory.usage_in_bytes",
                )
            } else {
                continue;
            };

            let root = Path::new(root);
            let mut dir: PathBuf = root.join(path.trim_start_matches('/'));
            loop {
                // a v2 limit of "max" does not parse, which
                // is the same as having no limit here
                if let (Some(limit), Some(usage)) = (
                    read_u64(&dir.join(limit_file)),
                    read_u64(&dir.join(usage_file)),
                ) {
                    let gap = limit.saturating_sub(usage);
                    headroom = Some(headroom.map_or(gap, |h| h.min(gap)));
                }
                if dir == root || !dir.pop() {
                    break;
                }
            }
        }
        headroom
    }

    fn read_u64(path: &Path) -> Option<u64> {
        fs::read_to_string(path).ok()?.trim().parse().ok()
    }
}

#[cfg(test)]
mod tests {
    #[test]
    #[cfg(target_os = "linux")]
    fn test_available_memory() {
        let available = super::available_memory().expect("known on linux");
        assert!(available > 0);
    }
}
