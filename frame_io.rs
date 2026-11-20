//! One-frame lookahead runs on the input worker; no buffers scale with video length.
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Backend {
    Auto,
    Threaded,
    Uring,
}

pub struct FrameReader {
    // Close the ring before the file, including on an I/O error.
    #[cfg(target_os = "linux")]
    ring: Option<io_uring::IoUring>,
    file: File,
    offset: u64,
    bytes: Option<Vec<u8>>,
}

impl FrameReader {
    pub fn new(file: File, offset: u64, size: usize, backend: Backend) -> io::Result<Self> {
        #[cfg(target_os = "linux")]
        let ring = if backend == Backend::Threaded {
            None
        } else {
            let setup = (|| {
                let ring = io_uring::IoUring::new(2)?;
                let mut probe = io_uring::Probe::new();
                ring.submitter().register_probe(&mut probe)?;
                if !probe.is_supported(io_uring::opcode::Read::CODE) {
                    return Err(io::Error::new(
                        io::ErrorKind::Unsupported,
                        "io_uring READ is unavailable",
                    ));
                }
                Ok(ring)
            })();
            match setup {
                Ok(ring) => Some(ring),
                Err(error) if backend == Backend::Auto => {
                    eprintln!("io_uring unavailable ({error}); using threaded input");
                    None
                }
                Err(error) => return Err(error),
            }
        };
        #[cfg(not(target_os = "linux"))]
        if backend == Backend::Uring {
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "io_uring requires Linux",
            ));
        }
        Ok(Self {
            #[cfg(target_os = "linux")]
            ring,
            file,
            offset,
            bytes: Some(vec![0; size]),
        })
    }

    pub fn name(&self) -> &'static str {
        #[cfg(target_os = "linux")]
        if self.ring.is_some() {
            return "io_uring";
        }
        "threaded"
    }

    pub fn read_frame(&mut self) -> io::Result<&[u8]> {
        #[cfg(target_os = "linux")]
        if let Some(ring) = &mut self.ring {
            use io_uring::{opcode, types};
            use std::os::fd::AsRawFd;
            let bytes = self
                .bytes
                .take()
                .ok_or_else(|| io::Error::other("reader failed previously"))?;
            let mut buffer = PendingBuffer {
                bytes: Some(bytes),
                in_flight: false,
            };
            let mut done = 0;
            while done < buffer.bytes.as_ref().unwrap().len() {
                let bytes = buffer.bytes.as_mut().unwrap();
                // READ completions use signed i32 byte counts; chunk exceptionally large frames.
                let count = (bytes.len() - done).min(i32::MAX as usize) as u32;
                let entry = opcode::Read::new(
                    types::Fd(self.file.as_raw_fd()),
                    bytes[done..].as_mut_ptr(),
                    count,
                )
                .offset(self.offset + done as u64)
                .build();
                // SAFETY: the Vec allocation remains owned by PendingBuffer and is never
                // accessed or resized until its completion is consumed. There is at most
                // one pending operation. On an unrecoverable wait error the guard keeps
                // the allocation alive, even while the ring cancels the outstanding read.
                unsafe { ring.submission().push(&entry) }
                    .map_err(|_| io::Error::other("io_uring submission queue unexpectedly full"))?;
                buffer.in_flight = true;
                let result = loop {
                    match ring.submit_and_wait(1) {
                        Ok(_) => {}
                        Err(error) if error.kind() == io::ErrorKind::Interrupted => continue,
                        Err(error) => return Err(error),
                    }
                    if let Some(completion) = ring.completion().next() {
                        break completion.result();
                    }
                };
                buffer.in_flight = false;
                if result < 0 {
                    let error = io::Error::from_raw_os_error(-result);
                    if error.kind() == io::ErrorKind::Interrupted {
                        continue;
                    }
                    return Err(error);
                }
                if result == 0 {
                    return Err(io::Error::new(
                        io::ErrorKind::UnexpectedEof,
                        "simulation was truncated during rendering",
                    ));
                }
                done += result as usize;
            }
            self.offset += done as u64;
            self.bytes = buffer.bytes.take();
            return Ok(self.bytes.as_ref().unwrap());
        }
        let bytes = self.bytes.as_mut().unwrap();
        self.file.seek(SeekFrom::Start(self.offset))?;
        self.file.read_exact(bytes)?;
        self.offset += bytes.len() as u64;
        Ok(bytes)
    }
}

#[cfg(target_os = "linux")]
struct PendingBuffer {
    bytes: Option<Vec<u8>>,
    in_flight: bool,
}

#[cfg(target_os = "linux")]
impl Drop for PendingBuffer {
    fn drop(&mut self) {
        if self.in_flight {
            // A fatal submission/wait failure does not establish that the kernel has
            // stopped using the address. Leak this one buffer rather than allow UAF.
            // The pipeline immediately aborts; ordinary I/O errors consume their CQE
            // and free the buffer normally.
            if let Some(bytes) = self.bytes.take() {
                std::mem::forget(bytes);
            }
        }
    }
}
