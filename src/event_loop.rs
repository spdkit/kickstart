// imports

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*imports][imports:1]]
use std::path::PathBuf;

use crate::common::*;
use crate::core::*;
// imports:1 ends here

// ctrlc event
// When user presses Ctrl-C, send signal using signal channel.

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*ctrlc%20event][ctrlc event:1]]
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

// runcmd event
// 搜索流程要做的事
// 1. 执行搜索流程
// 2. 响应用户控制信号.

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*runcmd%20event][runcmd event:1]]
type Job = Option<String>;

use std::sync::atomic::Ordering;
fn runcmd_channel(
    control: std::sync::Arc<JobFlag>,
) -> Result<(cbchan::Sender<Job>, cbchan::Receiver<Job>)> {
    let (sender, receiver) = cbchan::unbounded();

    let t = std::time::Duration::from_secs(2);
    let _sender = sender.clone();
    std::thread::spawn(move || loop {
        match JobType::from(&control) {
            JobType::Run => {
                let seeds = crate::search::prepare_seeds();
                let _ctrl = control.clone();
                crate::search::evolve_from_seeds(seeds, &_ctrl).unwrap();
            }
            JobType::Edit => {
                let msg = format!("i am ready!");
                let _ = _sender.send(Some(msg)).unwrap();
                std::thread::sleep(t);
            }
            JobType::Stop => {
                let _ = _sender.send(None).unwrap();
                break;
            }
            _ => unimplemented!(),
        }
    });

    Ok((sender, receiver))
}
// runcmd event:1 ends here

// pub/run
// 1. 在后台线程中执行搜索流程.
// 2. 监视用户ctrl-c信号, 收到后进入REPL
// 3. 在REPL界面, 用户发控制信号给搜索流程.

// [[file:~/Workspace/Programming/structure-predication/kickstart/kickstart.note::*pub/run][pub/run:1]]
pub fn enter_main_loop() -> Result<()> {
    let pid = std::process::id();
    println!("Issue cmd `kill -USR1 {}` to interrupt ...", pid);

    let flag = JobFlag::new(JobType::Run.flag());
    let control_flag = std::sync::Arc::new(flag);
    let signal_events = signal_channel(control_flag.clone())?;
    let (cmd_tx, cmd_rx) = runcmd_channel(control_flag.clone())?;

    loop {
        cbchan::select! {
            // wait for user interruption
            recv(signal_events) -> sig => {
                match sig {
                    Ok(10) | Ok(SIGINT) => {
                        let stop = crate::interactive::start_repl(control_flag.clone())?;
                        if stop {
                            info!("Terminating now ...");
                            break;
                        } else {
                            info!("Continue running");
                        }
                    }
                    Ok(sig) => {
                        warn!("unprocessed signal {:?}", sig);
                    }
                    Err(e) => {
                        error!("Process signal hook failed: {:?}", e);
                        break;
                    }
                }
            }
            recv(cmd_rx) -> msg => {
                match msg {
                    Ok(Some(Job)) => {
                        info!("find new job");
                    }
                    Ok(None) => {
                        info!("Job completed.");
                        break;
                    }
                    Err(e) => {
                        error!("found error when running cmd:\n {}", e);
                        break;
                    }
                }
            }
        }
    }

    Ok(())
}
// pub/run:1 ends here
