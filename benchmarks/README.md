# Benchmarks

Performance measurement scripts for the QBastion Falcon signing pipeline.

## Prerequisites

```bash
cd hydra/falcon_py
pip install pycryptodome numpy beartype
```

## Scripts

### `bench_falcon_signing.py` — End-to-end sign/verify latency

Measures `sign()` and `verify()` latency across all Falcon parameter sets
(n = 64, 128, 256, 512, 1024).

```bash
python benchmarks/bench_falcon_signing.py
```

Example output:
```
QBastion — Falcon Signing Benchmark
------------------------------------------------------------

Falcon-512  (n=512, sig_bytelen=666)
  sign()                  22.510 ms  (20 iterations)
  verify()                 0.183 ms  (20 iterations)

Falcon-1024  (n=1024, sig_bytelen=1280)
  sign()                 105.381 ms  (10 iterations)
  verify()                 0.381 ms  (10 iterations)
```

### `bench_sampler.py` — Gaussian sampler throughput

Profiles the `samplerz` discrete Gaussian sampler, which is the dominant
cost in Falcon-1024 signing.

```bash
python benchmarks/bench_sampler.py
```

Example output:
```
QBastion — SamplerZ Throughput Benchmark
------------------------------------------------------------
       sigma     mean (µs)    stdev (µs)     samples/s
------------------------------------------------------------
      1.3000         2.318         0.503       431,469
      1.5000         2.287         0.491       437,246
      1.7000         2.265         0.477       441,530
      1.8204         2.301         0.482       434,620
```

## Notes

- All timings are wall-clock (not CPU time) and include Python interpreter overhead.
- For production throughput estimations, the Cairo-based verifier benchmarks
  are more relevant — see `make falcon-execute` and `make falcon-prove` in the root Makefile.
