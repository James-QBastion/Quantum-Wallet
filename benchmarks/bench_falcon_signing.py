"""
Benchmarking script for Falcon signature generation throughput.

Measures sign() and verify() latency across all supported Falcon parameter
sets (n = 64 ... 1024). Results are printed in a table grouped by degree.

Usage:
    cd hydra/falcon_py
    python ../../benchmarks/bench_falcon_signing.py
"""
import sys
import os
import time

# Allow importing from hydra/falcon_py when running from the project root
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "hydra", "falcon_py"))

from falcon import Falcon, params  # noqa: E402

# ── Configuration ─────────────────────────────────────────────────────────────

# Message used in all benchmarks (realistic transaction payload size)
BENCH_MESSAGE = b"qbastion:bench:" + b"A" * 64

# Number of sign/verify iterations per parameter set
ITERATIONS = {
    64:   200,
    128:  100,
    256:  50,
    512:  20,
    1024: 10,
}

SEPARATOR = "-" * 60


# ── Helpers ───────────────────────────────────────────────────────────────────

def bench(label: str, fn, iters: int) -> float:
    """Run `fn` `iters` times and return the average milliseconds per call."""
    start = time.perf_counter()
    for _ in range(iters):
        fn()
    elapsed = time.perf_counter() - start
    ms_per_call = (elapsed / iters) * 1000
    print(f"  {label:<20} {ms_per_call:>10.3f} ms  ({iters} iterations)")
    return ms_per_call


# ── Main benchmark loop ───────────────────────────────────────────────────────

def main():
    print("QBastion — Falcon Signing Benchmark")
    print(SEPARATOR)

    for degree in sorted(params.keys()):
        if degree < 64:
            continue  # Skip tiny parameter sets — not used in practice

        iters = ITERATIONS.get(degree, 10)
        print(f"\nFalcon-{degree}  (n={degree}, sig_bytelen={params[degree].sig_bytelen})")

        prover = Falcon(degree)

        # Key generation (done once — not part of hot path)
        sk, vk = prover.keygen()

        # Sign benchmark
        bench(
            "sign()",
            lambda: prover.sign(sk, BENCH_MESSAGE),
            iters,
        )

        # Generate one valid signature for the verify benchmark
        sig = prover.sign(sk, BENCH_MESSAGE)

        # Verify benchmark
        bench(
            "verify()",
            lambda: prover.verify(vk, BENCH_MESSAGE, sig),
            iters,
        )

    print(f"\n{SEPARATOR}")
    print("Done.")


if __name__ == "__main__":
    main()
