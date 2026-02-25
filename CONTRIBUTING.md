# Contributing to QBastion

Thank you for your interest in contributing to the QBastion quantum-resistant node protection project! This document provides everything you need to get started.

---

## Table of Contents

- [Development Setup](#development-setup)
- [Project Structure](#project-structure)
- [Branching Strategy](#branching-strategy)
- [Commit Convention](#commit-convention)
- [Code Quality](#code-quality)
- [Submitting a Pull Request](#submitting-a-pull-request)

---

## Development Setup

### Cairo / Scarb

This project uses [Scarb](https://docs.swmansion.com/scarb/) for Cairo package management.

```bash
# Install via asdf (recommended — versions are pinned in .tool-versions)
asdf install

# Verify
scarb --version
```

### Python (hydra module)

The `hydra/` directory contains Python reference implementations.

```bash
cd hydra/falcon_py
pip install pycryptodome numpy beartype pre-commit

# Run tests
python test.py
```

### Pre-commit hooks

```bash
pip install pre-commit
pre-commit install
```

---

## Project Structure

```
packages/falcon/     # Cairo Falcon-512 verification package
hydra/               # Python reference implementations
  falcon_py/         # Falcon signing / KAT generator
  bounded_int_circuit/ # Cairo code generator
tests/               # Python test suite
docs/plans/          # Design notes and implementation plans
```

---

## Branching Strategy

| Branch pattern | Purpose |
|----------------|---------|
| `main` | Stable, always passing CI |
| `feat/<name>` | New feature |
| `fix/<name>` | Bug fix |
| `docs/<name>` | Documentation only |
| `refactor/<name>` | Code refactoring, no functional change |

Always branch off `main` and open a PR back to `main`.

---

## Commit Convention

We follow [Conventional Commits](https://www.conventionalcommits.org/):

```
<type>(<scope>): <short description>

[optional body]
[optional footer]
```

**Types:** `feat`, `fix`, `docs`, `refactor`, `test`, `chore`, `perf`

**Scopes:** `cairo`, `falcon-py`, `ci`, `docs`, `ntru`, `ntt`

**Examples:**

```
feat(cairo): add Sphincs+ batch verifier entrypoint
fix(falcon-py): correct sigma_max boundary in samplerz
docs: update roadmap with Stwo proving benchmarks
```

---

## Code Quality

### Cairo

```bash
scarb fmt          # autoformat
scarb fmt --check  # lint only
```

### Python

Pre-commit handles `black`, `isort`, and `flake8` automatically on every commit. 
To run manually:

```bash
pre-commit run --all-files
```

---

## Submitting a Pull Request

1. Fork the repository and create a feature branch.
2. Make your changes, including tests where relevant.
3. Ensure CI passes locally (`scarb test`, `python test.py`).
4. Open a PR with a clear description of what changed and why.
5. Address any review comments promptly.

We aim to review PRs within 3 business days.
