# Changelog

All notable changes to QBastion will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [Unreleased]

### Added
- `docs/architecture.md` — high-level QBastion system design document
- `benchmarks/` — profiling scripts for signature generation throughput
- `CHANGELOG.md` — this file
- `CONTRIBUTING.md` — developer onboarding guide
- `SECURITY.md` — vulnerability disclosure policy
- `.pre-commit-config.yaml` — automated code quality checks (black, isort, flake8)
- `.github/workflows/cairo-ci.yml` — Cairo build and test CI pipeline
- `.github/workflows/python-tests.yml` — Python Falcon test suite CI pipeline
- Type annotations added to `ntrugen.py`
- Comprehensive docstrings and algorithm references across  
  `ffsampling.py`, `samplerz.py`, `fft.py`, `ntt.py`, `common.py`, `encoding.py`, `rng.py`

### Changed
- `hydra/falcon_py/` refactored with private helper names and QBastion-specific module docstrings
- `hydra/falcon_py/README.md` rewritten with table-based module overview

---

## [0.1.0] — 2025-01-01

### Added
- Falcon-512 signature verification circuit in Cairo (`packages/falcon`)
- Sphincs+-128s with SHA-2 (simple mode) verification
- Sphincs+-128s with BLAKE2s and 4-byte aligned address encoding
- Python reference implementation of Falcon (`hydra/falcon_py/`)
- NTT-based polynomial multiplication in Cairo (`packages/falcon/src/ntt_constants.cairo`)
- Ingress gateway module skeleton (`packages/falcon/src/ingress.cairo`)
- BoundedInt circuit compiler (`hydra/`)
- NTRU key generation and polynomial arithmetic (`hydra/falcon_py/ntrugen.py`)
- Prover parameter configuration (`prover_params.json`)

### References
- [BIP360](https://bip360.org/) — Pay to Quantum Resistant Hash
- [PQ Signatures and Scaling Bitcoin with STARKs](https://delvingbitcoin.org/t/post-quantum-signatures-and-scaling-bitcoin-with-starks/1584)
