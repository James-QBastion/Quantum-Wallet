# QBastion Falcon Python Integration

A Python reference implementation of the Falcon post-quantum lattice signature scheme,
adapted for the QBastion node protection pipeline. This module is used mainly for:

- Generating and verifying Falcon-512 / Falcon-1024 signatures during testing
- Producing KAT (Known Answer Test) vectors for cross-validation against the Cairo verifier
- Benchmarking signature generation throughput on commodity hardware

The underlying algorithm is the Falcon lattice-based signature scheme
(https://falcon-sign.info/), specifically NIST security levels I and V.

---

## Structure

| File | Purpose |
|------|---------|
| `common.py` | Modular constants, `poly_split` / `poly_merge` helpers |
| `samplerz.py` | Discrete Gaussian sampler over Z (RCDT + Bernoulli exp) |
| `fft.py` | Anti-cyclotomic FFT/IFFT over R[x]/(x^n+1) |
| `ntt.py` | NTT/INTT over Z_q[x]/(x^n+1), q=12289 |
| `rng.py` | ChaCha20 deterministic byte generator |
| `ffsampling.py` | Fast Fourier nearest-plane (ffNP) sampler |
| `ntrugen.py` | NTRU key generation |
| `encoding.py` | Signature compression/decompression codec |
| `falcon.py` | Top-level `Falcon` class: keygen, sign, verify |
| `test.py` | End-to-end test suite with KAT validation |

---

## Quick Start

```python
from falcon import Falcon

prover = Falcon(512)             # Falcon-512
sk, vk = prover.keygen()
sig    = prover.sign(sk, b"my message")
assert prover.verify(vk, b"my message", sig)
```

Keys and signatures are byte strings, making them easy to serialize and pass
to the Cairo verifier.

---

## Running Tests

```bash
make test
```

Expected output (approximate timings vary by machine):

```
Test Sig KATs       : OK
Test SamplerZ KATs  : OK

Test battery for n = 64
Test FFT            : OK
Test NTT            : OK
Test NTRUGen        : OK
Test ffNP           : OK
Test Compress       : OK
Test Signature      : OK
```

---

## Profiling

Run `make profile` to generate a call-graph profile using `pyprof2calltree`
and `kcachegrind`.

---

## Author

**James Brown**

---

## License

MIT — see `LICENSE`.
