<!--
Thanks for contributing to phylo-rs. See CONTRIBUTING.md for build, test, and
review details.
-->

## What does this change?

<!-- A short description of the change and why it is needed. Link any related
issue with "Fixes #123". -->

## How was it tested?

<!-- Which tests did you add or run? For a bug fix, ideally a test that fails
before this change. -->

## Checklist

- [ ] `cargo fmt --all` has been run
- [ ] `cargo clippy --all-targets --all-features` is clean
- [ ] `cargo test` passes
- [ ] Tests added or updated for the change
- [ ] Touched parallel code paths? Also ran `cargo test --features parallel`
- [ ] Public API changed? Doc comments updated

<!--
CI runs fmt, clippy, build, test, and bench --no-run, and must pass before
merge. Your branch also needs to be up to date with main.
-->
