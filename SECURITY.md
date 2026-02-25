# Security Policy

## Supported Versions

This project is currently in active development. Security fixes are applied to the `main` branch only.

| Version | Supported |
|---------|-----------|
| `main` | ✅ |
| Historical tags | ❌ |

---

## Scope

This security policy applies to:

- The **Cairo verifier** packages under `packages/`
- The **Python reference implementation** under `hydra/falcon_py/`
- Any **cryptographic primitives** (Falcon, Sphincs+, NTRU routines)

**Out of scope:**

- Known limitations of the Falcon specification itself
- Issues in third-party dependencies (report those upstream)
- Build tooling configuration issues that have no security impact

---

## Reporting a Vulnerability

> ⚠️ **Please do NOT open a public GitHub issue for security vulnerabilities.**

To report a security issue privately:

1. Email **security@george-quantum.dev** with the subject line: `[SECURITY] <brief description>`
2. Include:
   - A clear description of the vulnerability
   - Steps to reproduce or a proof-of-concept
   - The affected file(s) and line numbers if known
   - Your assessment of severity (critical / high / medium / low)

We will acknowledge your report within **48 hours** and provide a fix timeline within **7 days**.

---

## Cryptographic Notes

This project implements post-quantum cryptographic primitives. Note that:

- The Python code in `hydra/falcon_py/` is a **reference implementation** — it is **not** hardened for production use or side-channel resistance.
- The Cairo implementation is designed for **STARK-based batch verification** in a zkVM context, and its security depends on the underlying Stwo prover's soundness.

---

## Acknowledgements

We thank security researchers who responsibly disclose vulnerabilities. Confirmed reporters will be credited in release notes (unless they prefer anonymity).
