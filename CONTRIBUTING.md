# Contributing to phylo-rs

Thanks for your interest in contributing. This document covers how to build,
test, and submit changes. For what the library does and how to use it, see the
[README](README.md).

By participating you agree to abide by our [Code of Conduct](CODE_OF_CONDUCT.md).

## Getting started

```bash
git clone https://github.com/sriram98v/phylo-rs
cd phylo-rs
cargo build
cargo test
```

The minimum supported Rust version is **1.87**, declared as `rust-version` in
`Cargo.toml` and checked by CI against a freshly resolved dependency tree.
Day-to-day development happens on stable.

Note that 1.87 is imposed by dependencies, not by this crate: the library
itself compiles on 1.80. `vers-vecs` declares no `rust-version` of its own, so
cargo cannot filter it by MSRV and it resolves to the newest `1.x`, which
currently needs 1.87. If that dependency raises its requirements again, the
`Check MSRV` CI job will fail and `rust-version` needs to move with it.

## Development commands

```bash
cargo build                                # build the library
cargo test                                 # run all tests
cargo test --test tree-tests <test_name>   # run a single integration test
cargo fmt --all                            # format (required, see below)
cargo clippy --all-targets --all-features  # lint
cargo bench --no-run                       # build benchmarks without running
cargo bench                                # run benchmarks (divan)
cargo run --example phylogenetic-diversity # run an example
```

## Feature flags

| Flag | Default | Effect |
| --- | --- | --- |
| `non_crypto_hash` | yes | use `fxhash` instead of `std` hashing |
| `simple_rooted_tree` | yes | the `SimpleRootedTree` implementation and its aliases |
| `parallel` | no | rayon-based parallel computation |
| `serde` | no | serialization support |

If your change touches parallel code paths, test both ways:

```bash
cargo test
cargo test --features parallel
```

### Hashing is feature-gated

This trips people up. When `non_crypto_hash` is enabled (the default),
`FxHashMap`/`FxHashSet` stand in for the `std` types. Code that uses hash maps
must gate the import rather than importing from `std` directly:

```rust
#[cfg(feature = "non_crypto_hash")]
use fxhash::{FxHashMap as HashMap, FxHashSet as HashSet};
#[cfg(not(feature = "non_crypto_hash"))]
use std::collections::{HashMap, HashSet};
```

Importing `std::collections::HashMap` unconditionally will compile, but it
diverges from the rest of the codebase and defeats the flag.

## How the code is organised

The design is trait-based: narrow traits compose to provide tree functionality,
and the crate ships one concrete implementation plus type aliases (`PhyloTree`,
`PhyloNode`). `use phylo::prelude::*;` is the standard entry point and imports
all public traits and types.

**Nodes** (`src/node/`) — `Node<T, W, Z>` is generic over taxa type, edge
weight, and node annotation. The traits build up incrementally:

- `RootedTreeNode` — id, parent, children
- `RootedMetaNode` — taxa annotations
- `RootedWeightedNode` — edge weights
- `RootedZetaNode` — numeric annotations, used by distance metrics

**Trees** (`src/tree/`) — `SimpleRootedTree<T, W, Z>` is arena-allocated
(`Vec<Option<Node<T,W,Z>>>`) and indexed by `NodeID`:

- `simple_rtree.rs` — the `RootedTree` / `RootedMetaTree` / `RootedWeightedTree`
  traits and the `SimpleRootedTree` struct
- `ops.rs` — mutations: SPR, NNI, reroot, subtree extraction, contraction
- `distances.rs` — RF, weighted RF, cluster affinity, cophenetic, distance matrices
- `io.rs` — Newick and Nexus parsing and serialization
- `simulation.rs` — random tree generation (Yule, uniform)

**Iteration** (`src/iter/node_iter.rs`) — `DFS`, `BFS`, `PreOrder`, `Ancestors`,
`EulerWalk`, `Clusters`, implemented on any `RootedTree`. `EulerWalk` supports
O(1) LCA queries via a borrowing `LcaOracle` built with `tree.lca()`.

New code should land in the module that owns that concern, and should be generic
over the traits rather than reaching for `SimpleRootedTree` directly where a
trait bound would do.

## Tests

Integration tests live in `tests/tree-tests.rs` and `tests/node-tests.rs`.
Benchmarks use [divan](https://docs.rs/divan) in `benches/main.rs`.

Please include a test with any bug fix or new feature. A bug fix without a test
that fails before the change is hard to review and easy to regress.

## Before you open a pull request

CI runs exactly these, and all of them must pass:

```bash
cargo fmt --all --check
cargo clippy --all-targets --all-features
cargo build
cargo test
cargo bench --no-run
```

Formatting is enforced, not advisory — run `cargo fmt --all` before you commit.

## Pull requests

`main` is protected. Every change goes through a pull request that:

- passes CI (the `Lint, build, test, and build benchmarks` job is required)
- is up to date with `main` — rebase if `main` has moved
- has an approving review from a code owner (see [CODEOWNERS](.github/CODEOWNERS))

Keep unrelated changes in separate pull requests. In particular, avoid mixing
reformatting with behavioural changes — it makes review much harder.

### Commit messages

We follow [Conventional Commits](https://www.conventionalcommits.org): a `type:`
prefix and a short imperative summary.

```
fix: correct branch length parsing for zero-length edges

Longer explanation of the why, wrapped at 72 characters, if the summary
is not self-explanatory.
```

Types in use: `feat`, `fix`, `refactor`, `docs`, `test`, `chore`, `perf`, `ci`.
This is a convention rather than a gate — nothing in CI validates it — but
following it keeps the history readable.

## Reporting bugs and requesting features

Open an [issue](https://github.com/sriram98v/phylo-rs/issues). For bugs, a
Newick string or short snippet that reproduces the problem is worth more than a
description of it. Note that issues are closed automatically after 30 days of
inactivity.

## Licence

phylo-rs is MIT licensed. By contributing, you agree that your contributions
will be licensed under the same terms. See [LICENSE](LICENSE).

## Citation

If you use phylo-rs in academic work, please cite the paper linked in the
[README](README.md#citation).
