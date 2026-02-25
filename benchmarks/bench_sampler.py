"""
Benchmarking script for the Gaussian sampler (samplerz) used in Falcon signing.

Measures throughput of the samplerz routine across a range of sigma values.
This is relevant because Gaussian sampling is the dominant cost in Falcon-1024
signing (~80% of total signing time for large parameter sets).

Usage:
    cd hydra/falcon_py
    python ../../benchmarks/bench_sampler.py
"""
import sys
import os
import time
from statistics import mean, stdev

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "hydra", "falcon_py"))

from samplerz import samplerz, SIGMA_MAX  # noqa: E402

SEPARATOR = "-" * 60
ITERATIONS = 50_000


def bench_sigma(sigma: float, iters: int) -> tuple[float, float]:
    """Benchmark samplerz at a fixed sigma. Returns (mean_us, stdev_us)."""
    sigmin = 1.2
    mu     = 0.0
    times  = []
    for _ in range(iters):
        t0 = time.perf_counter()
        samplerz(mu, sigma, sigmin)
        times.append((time.perf_counter() - t0) * 1e6)
    return mean(times), stdev(times)


def main():
    print("QBastion — SamplerZ Throughput Benchmark")
    print(SEPARATOR)
    print(f"  {'sigma':>10}  {'mean (µs)':>12}  {'stdev (µs)':>12}  {'samples/s':>12}")
    print(SEPARATOR)

    test_sigmas = [1.3, 1.5, 1.7, SIGMA_MAX - 0.01]
    for sigma in test_sigmas:
        avg_us, sd_us = bench_sigma(sigma, ITERATIONS)
        samples_per_sec = 1_000_000 / avg_us
        print(f"  {sigma:>10.4f}  {avg_us:>12.3f}  {sd_us:>12.3f}  {samples_per_sec:>12,.0f}")

    print(SEPARATOR)
    print("Done.")


if __name__ == "__main__":
    main()
