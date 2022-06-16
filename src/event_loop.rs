// [[file:../kickstart.note::*imports][imports:1]]
use std::path::PathBuf;

use crate::common::*;
use crate::core::*;
// imports:1 ends here

// [[file:../kickstart.note::*ctrlc event][ctrlc event:1]]
use crossbeam_channel as cbchan;
use signal_hook::{SIGINT, SIGUSR1};

type Result<T> = std::result::Result<T, Error>;

fn signal_channel(control_flag: std::sync::Arc<JobFlag>) -> Result<cbchan::Receiver<i32>> {
    let signals = signal_hook::iterator::Signals::new(&[SIGINT, SIGUSR1])?;
    let (sender, receiver) = cbchan::bounded(1);

    std::thread::spawn(move || {
        // It will wait for signal forever
        for sig in &signals {
            let _ = sender.send(sig);
        }
    });

    Ok(receiver)
}
// ctrlc event:1 ends here
