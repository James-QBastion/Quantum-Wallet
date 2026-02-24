"""
Discrete Gaussian sampler for the QBastion lattice signer.

Implements an integer Gaussian sampler D_{Z, mu, sigma} using:
  1. A half-Gaussian base sampler via a reverse cumulative distribution table (RCDT)
  2. A Bernoulli exponential test to accept/reject samples

The RCDT table and exp-polynomial coefficients are derived from the Falcon
specification and the FACCT paper (https://doi.org/10.1109/TC.2019.2940949).
"""
from math import floor
from os import urandom as _system_random


# ── Gaussian parameters ───────────────────────────────────────────────────────

# Hard upper bound on the per-coordinate sigma across all Falcon parameter sets
SIGMA_MAX = 1.8205

# Precomputed reciprocal: 1 / (2 * sigma_max^2), used in rejection test
_HALF_INV_SIGMA2 = 1.0 / (2.0 * SIGMA_MAX ** 2)

# Alias kept for backward-compat imports
MAX_SIGMA    = SIGMA_MAX
INV_2SIGMA2  = _HALF_INV_SIGMA2

# Bit-precision of the RCDT random draw
_RCDT_BITS = 72
RCDT_PREC  = _RCDT_BITS  # legacy alias

# Natural log of 2 and its reciprocal
_LN2  = 0.69314718056
_ILN2 = 1.44269504089
LN2   = _LN2   # legacy alias
ILN2  = _ILN2  # legacy alias


# ── Reverse Cumulative Distribution Table ─────────────────────────────────────
# Precomputed threshold table for sampling a half-Gaussian with std dev SIGMA_MAX.
# Entry i represents the probability mass for z0 >= i (scaled to 2^72).
_GAUSS_RCDT = [
    3024686241123004913666,
    1564742784480091954050,
    636254429462080897535,
    199560484645026482916,
    47667343854657281903,
    8595902006365044063,
    1163297957344668388,
    117656387352093658,
    8867391802663976,
    496969357462633,
    20680885154299,
    638331848991,
    14602316184,
    247426747,
    3104126,
    28824,
    198,
    1,
]
RCDT = _GAUSS_RCDT  # legacy alias


# ── Polynomial approximation coefficients for exp(-x) ────────────────────────
# Lifted from FACCT (https://doi.org/10.1109/TC.2019.2940949).
# Approximates 2^-63 * sum(C[12 - i] * x^i) ≈ exp(-x) for 0 <= x < ln2.
_EXP_POLY_COEFFS = [
    0x00000004741183A3,
    0x00000036548CFC06,
    0x0000024FDCBF140A,
    0x0000171D939DE045,
    0x0000D00CF58F6F84,
    0x000680681CF796E3,
    0x002D82D8305B0FEA,
    0x011111110E066FD0,
    0x0555555555070F00,
    0x155555555581FF00,
    0x400000000002B400,
    0x7FFFFFFFFFFF4800,
    0x8000000000000000,
]
C = _EXP_POLY_COEFFS  # legacy alias


# ── Core sampling functions ───────────────────────────────────────────────────

def _sample_half_gaussian(rng=_system_random):
    """Draw z0 from {0, 1, ..., 18} according to 1/2 * D_{Z+, 0, SIGMA_MAX}.

    Uses the RCDT to implement a constant-time comparison against each
    cumulative threshold.

    Args:
        rng: callable that returns k random bytes (default: os.urandom)

    Returns:
        An integer z0 in [0, 18].
    """
    rand_int = int.from_bytes(rng(_RCDT_BITS >> 3), "little")
    z0 = 0
    for threshold in _GAUSS_RCDT:
        z0 += int(rand_int < threshold)
    return z0


# Legacy alias
basesampler = _sample_half_gaussian


def _approx_scaled_exp(x, ccs):
    """Compute an integer approximation of 2^63 * ccs * exp(-x).

    Uses a degree-12 polynomial approximation of exp(-x) built from
    _EXP_POLY_COEFFS. Both inputs must be strictly positive.

    Args:
        x:   positive floating-point exponent
        ccs: positive scaling factor

    Returns:
        Integer approximation of 2^63 * ccs * exp(-x).
    """
    acc = _EXP_POLY_COEFFS[0]
    x_fixed = int(x * (1 << 63))  # convert to fixed-point Q1.63
    for coeff in _EXP_POLY_COEFFS[1:]:
        acc = coeff - ((x_fixed * acc) >> 63)
    scale = int(ccs * (1 << 63)) << 1
    return (scale * acc) >> 63


# Legacy alias
approxexp = _approx_scaled_exp


def _bernoulli_exp_test(x, ccs, rng=_system_random):
    """Bernoulli trial: return 1 with probability ≈ ccs * exp(-x).

    Reads random bytes byte-by-byte and compares them MSB-first against the
    approximated probability, stopping as soon as the outcome is determined.
    Both inputs must be positive.

    Args:
        x:   positive floating-point exponent value
        ccs: positive scaling factor
        rng: source of random bytes

    Returns:
        True (accept) or False (reject)
    """
    # Decompose x = s * ln2 + r so that exp(-x) = 2^-s * exp(-r)
    shift = int(x * _ILN2)
    remainder = x - shift * _LN2
    shift = min(shift, 63)
    prob = (_approx_scaled_exp(remainder, ccs) - 1) >> shift
    # Compare byte by byte from most-significant to least-significant
    for bit_offset in range(56, -8, -8):
        rand_byte = int.from_bytes(rng(1), "little")
        delta = rand_byte - ((prob >> bit_offset) & 0xFF)
        if delta:
            return delta < 0
    return False


# Legacy alias
berexp = _bernoulli_exp_test


def samplerz(mu, sigma, sigmin, randombytes=_system_random):
    """Sample an integer z from the discrete Gaussian D_{Z, mu, sigma}.

    Implements the rejection sampler from the Falcon specification:
    draws a candidate from the half-Gaussian via the RCDT, promotes it to a
    full Gaussian candidate, then accepts or rejects via a Bernoulli exp test.

    Constraints: 1 < sigmin < sigma < SIGMA_MAX.

    Args:
        mu:          center of the distribution (float)
        sigma:       desired standard deviation (float)
        sigmin:      lower bound on sigma used for scaling (float)
        randombytes: entropy source, default os.urandom

    Returns:
        An integer z sampled from D_{Z, mu, sigma}.
    """
    floor_mu = int(floor(mu))
    frac_mu  = mu - floor_mu

    dss = 1.0 / (2.0 * sigma * sigma)
    ccs = sigmin / sigma

    while True:
        # Draw from the half-Gaussian
        z0 = _sample_half_gaussian(rng=randombytes)
        # Symmetrize: flip sign with probability 1/2
        sign_bit = int.from_bytes(randombytes(1), "little") & 1
        z = sign_bit + (2 * sign_bit - 1) * z0
        # Rejection test against the target Gaussian
        exp_arg  = ((z - frac_mu) ** 2) * dss
        exp_arg -= (z0 ** 2) * _HALF_INV_SIGMA2
        if _bernoulli_exp_test(exp_arg, ccs, rng=randombytes):
            return z + floor_mu
