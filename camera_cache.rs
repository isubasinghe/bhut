//! Small, versioned sidecar for exact camera bounds. Input data is never modified.
use std::fs::{self, File, OpenOptions};
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

pub type Bounds = ([f64; 2], f64);

pub fn signature(file: &File) -> io::Result<String> {
    let meta = file.metadata()?;
    let mut stamp = format!(
        "{} {:?} {:?}",
        meta.len(),
        meta.modified()?,
        meta.created().ok()
    );
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        // ctime also invalidates same-size edits whose mtime was restored; inode
        // identity detects replacement with another file at the same path.
        stamp.push_str(&format!(
            " {} {} {} {}",
            meta.dev(),
            meta.ino(),
            meta.ctime(),
            meta.ctime_nsec()
        ));
    }
    Ok(stamp)
}

pub struct CameraCache {
    pub path: PathBuf,
    signature: String,
    key: Vec<u8>,
}

impl CameraCache {
    pub fn new(
        input: &Path,
        signature: &str,
        particles: usize,
        frames: usize,
        angles: [f64; 3],
    ) -> io::Result<Self> {
        let source = fs::canonicalize(input)?;
        let mut path = source.as_os_str().to_os_string();
        path.push(".camera-cache");
        // Bump the version if projection, quantiles, or camera padding changes.
        // Resolution, codec, and render style do not affect these square bounds.
        let key = format!("bhut-camera-v1\n{source:?}\n{signature}\n{particles} {frames}\n{:016x} {:016x} {:016x}\n",
            angles[0].to_bits(), angles[1].to_bits(), angles[2].to_bits()).into_bytes();
        Ok(Self {
            path: path.into(),
            signature: signature.into(),
            key,
        })
    }

    pub fn unchanged(&self, file: &File) -> bool {
        signature(file).is_ok_and(|current| current == self.signature)
    }

    pub fn load(&self) -> Option<Bounds> {
        let expected = self.key.len() + 32;
        let mut bytes = Vec::new();
        File::open(&self.path)
            .ok()?
            .take(expected as u64 + 1)
            .read_to_end(&mut bytes)
            .ok()?;
        if bytes.len() != expected || !bytes.starts_with(&self.key) {
            return None;
        }
        let payload = &bytes[self.key.len()..];
        let stored = u64::from_le_bytes(payload[24..].try_into().ok()?);
        if checksum(&bytes[..expected - 8]) != stored {
            return None;
        }
        let values: Vec<_> = payload[..24]
            .chunks_exact(8)
            .map(|b| f64::from_le_bytes(b.try_into().unwrap()))
            .collect();
        if !values.iter().all(|v| v.is_finite()) || values[2] < 1e-12 {
            return None;
        }
        Some(([values[0], values[1]], values[2]))
    }

    pub fn store(&self, bounds: Bounds) -> io::Result<()> {
        let mut bytes = self.key.clone();
        for value in [bounds.0[0], bounds.0[1], bounds.1] {
            bytes.extend(value.to_le_bytes());
        }
        bytes.extend(checksum(&bytes).to_le_bytes());
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let mut temp = self.path.as_os_str().to_os_string();
        temp.push(format!(
            ".tmp-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        let temp = PathBuf::from(temp);
        let mut file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temp)?;
        struct Cleanup(PathBuf);
        impl Drop for Cleanup {
            fn drop(&mut self) {
                let _ = fs::remove_file(&self.0);
            }
        }
        let _cleanup = Cleanup(temp.clone());
        file.write_all(&bytes)?;
        file.sync_all()?;
        drop(file);
        // Readers see the previous complete record or the new one, never a partial write.
        fs::rename(temp, &self.path)
    }
}

fn checksum(bytes: &[u8]) -> u64 {
    // Detect accidental cache corruption; this is not a security/content hash.
    bytes.iter().fold(0xcbf29ce484222325, |hash, &b| {
        (hash ^ u64::from(b)).wrapping_mul(0x100000001b3)
    })
}
