// imports

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*imports][imports:1]]
use std::path::PathBuf;

use crate::common::*;
// imports:1 ends here

// ctrl-c
// When user presses Ctrl-C, send signal using signal channel.

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*ctrl-c][ctrl-c:1]]
use crossbeam_channel as cbchan;
use signal_hook::SIGINT;

type Result<T> = std::result::Result<T, Error>;

fn signal_channel() -> Result<cbchan::Receiver<i32>> {
    let signals = signal_hook::iterator::Signals::new(&[SIGINT])?;
    let (sender, receiver) = cbchan::bounded(1);

    std::thread::spawn(move || {
        for sig in signals.forever() {
            let _ = sender.send(sig);
        }
    });

    Ok(receiver)
}
// ctrl-c:1 ends here

// cmd

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*cmd][cmd:1]]
fn runcmd_channel() -> Result<cbchan::Receiver<u32>> {
    let (sender, receiver) = cbchan::bounded(1);

    std::thread::spawn(move || {
        {
            // normal termination
            let _ = sender.send(0);
            unimplemented!();
        }
    });

    Ok(receiver)
}
// cmd:1 ends here

// run

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*run][run:1]]
fn enter_event_loop() -> Result<()> {
    let signal_events = signal_channel()?;
    let runcmd_events = runcmd_channel()?;

    println!("Press Ctrl-C to stop ...");
    let mut kill = false;
    loop {
        cbchan::select! {
            recv(signal_events) -> sig => {
                match sig {
                    Ok(SIGINT) => {
                        println!("Try to gracefully exit ...");
                        break;
                    }
                    Ok(sig) => {
                        warn!("unprocessed signal {:?}", sig);
                    }
                    Err(e) => {
                        eprintln!("Process signal hook failed: {:?}", e);
                        break;
                    }
                }
            }
            recv(runcmd_events) -> msg => {
                match msg {
                    Ok(0) => {
                        println!("Job completed.");
                        break;
                    }
                    Ok(pid) => {
                        info!("script session id: {}", pid);
                    }
                    _ => unreachable!()
                }
            }
        }
    }

    Ok(())
}
// run:1 ends here
