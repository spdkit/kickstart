# AGENTS.md

Operational notes for future coding agents working on this repository.

## Project positioning

- This repository is a public research artifact for SI.
- Distribution is binary-first (Linux x86_64 release assets).
- Source is available for method transparency; full local build parity is not required for SI users.

## Documentation policy

- Public-facing docs must be in English.
- Use Markdown only for maintained docs.
- `user-guide.md` is the canonical command guide.
- Do not add or maintain `.note`/Org-mode files for public docs.

## Release policy

- Release version is managed in `Cargo.toml`.
- Official binaries are `musl` static builds.
- Release binaries are uploaded as GitHub Release assets.
- Do **not** commit release binaries into git history.
- Attach checksums (`SHA256SUMS`) in each GitHub Release.

## Safety and hygiene

- Do not commit private paths, credentials, or internal infrastructure details.
- Keep `README.md` aligned with manuscript/SI wording.
- Before release, verify binaries run with:
  - `kickstart --help`
  - `kickgen --help`
  - `kickmut --help`
