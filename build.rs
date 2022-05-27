// [[file:kickstart.note::fb072ac9][fb072ac9]]
fn get_commit_hash() -> Option<String> {
    std::process::Command::new("git")
        .args(&["describe", "--always", "--dirty=-modified"])
        .output()
        .ok()
        .and_then(|r| String::from_utf8(r.stdout).ok())
}

fn get_commit_date() -> Option<String> {
    std::process::Command::new("git")
        .args(&["log", "-1", "--date=iso", "--pretty=format:%cd"])
        .output()
        .ok()
        .and_then(|r| String::from_utf8(r.stdout).ok())
}

fn main() {
    println!("cargo:rustc-env=GIT_HASH={}", get_commit_hash().unwrap_or_default());
    println!("cargo:rustc-env=COMMIT_DATE={}", get_commit_date().unwrap_or_default());
}
// fb072ac9 ends here
