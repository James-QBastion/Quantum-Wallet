"""
Fast Fourier Transform over R[x] / (x^n + 1) for the QBastion signing backend.

This module implements a recursive Cooley-Tukey FFT and its inverse (IFFT),
specifically for two-power-degree anti-cyclotomic rings. The transform is
needed to compute fast polynomial multiplication, division, and adjoint
operations used in the Falcon lattice signing procedure.

The split / merge helpers delegate to common.poly_split / poly_merge.
Pre-computed complex roots are loaded from fft_constants.roots_dict.
"""
from common import poly_split, poly_merge

# Import split/merge under legacy names so existing callers still work
split = poly_split
merge = poly_merge

from fft_constants import roots_dict  # complex n-th roots of unity, indexed by ring degree


# ── Radix-2 recursive butterflies (FFT domain) ───────────────────────────────

def _fft_split(poly_fft):
    """Split an FFT-represented polynomial into two half-degree FFT polynomials.

    Applies the splitfft_2 algorithm from the Falcon specification (Algorithm 1).

    Args:
        poly_fft: a polynomial in FFT representation (length n)

    Returns:
        [f0_fft, f1_fft] each of length n/2
    """
    deg = len(poly_fft)
    roots = roots_dict[deg]
    half = deg // 2
    f0_fft = [0.0] * half
    f1_fft = [0.0] * half
    for k in range(half):
        lo = poly_fft[2 * k]
        hi = poly_fft[2 * k + 1]
        f0_fft[k] = 0.5 * (lo + hi)
        f1_fft[k] = 0.5 * (lo - hi) * roots[2 * k].conjugate()
    return [f0_fft, f1_fft]


split_fft = _fft_split  # public alias


def _fft_merge(pair_fft):
    """Merge two half-degree FFT polynomials into a full-degree FFT polynomial.

    Applies the mergefft_2 algorithm from the Falcon specification (Algorithm 2).

    Args:
        pair_fft: [f0_fft, f1_fft], each of length n/2

    Returns:
        The reconstructed FFT polynomial of degree n.
    """
    f0_fft, f1_fft = pair_fft
    half = len(f0_fft)
    deg  = 2 * half
    roots = roots_dict[deg]
    result = [0.0] * deg
    for k in range(half):
        w_f1 = roots[2 * k] * f1_fft[k]
        result[2 * k]     = f0_fft[k] + w_f1
        result[2 * k + 1] = f0_fft[k] - w_f1
    return result


merge_fft = _fft_merge  # public alias


# ── Forward and inverse FFT ───────────────────────────────────────────────────

def fft(poly):
    """Compute the FFT of a polynomial in coefficient form.

    Recursively decomposes into degree-n/2 subproblems down to n == 2,
    where the base case is a single complex DFT butterfly.

    Args:
        poly: coefficient list of length n (must be a power of two, n >= 2)

    Returns:
        FFT representation (list of n complex values)
    """
    deg = len(poly)
    if deg > 2:
        even, odd       = poly_split(poly)
        even_fft        = fft(even)
        odd_fft         = fft(odd)
        return _fft_merge([even_fft, odd_fft])
    # Base case n == 2: two-point DFT
    return [poly[0] + 1j * poly[1],
            poly[0] - 1j * poly[1]]


def ifft(poly_fft):
    """Compute the inverse FFT, recovering a polynomial in coefficient form.

    Args:
        poly_fft: FFT representation (list of n complex values)

    Returns:
        Coefficient list of length n (values are real but returned as floats)
    """
    deg = len(poly_fft)
    if deg > 2:
        f0_fft, f1_fft = _fft_split(poly_fft)
        even           = ifft(f0_fft)
        odd            = ifft(f1_fft)
        return poly_merge([even, odd])
    # Base case n == 2
    return [poly_fft[0].real, poly_fft[0].imag]


# ── Polynomial arithmetic (coefficient representation) ───────────────────────

def add(f, g):
    """Pointwise addition of two polynomials (coefficient representation)."""
    assert len(f) == len(g)
    return [f[i] + g[i] for i in range(len(f))]


def neg(f):
    """Negation of a polynomial (works in any representation)."""
    return [-c for c in f]


def sub(f, g):
    """Subtraction of two polynomials (any representation)."""
    return add(f, neg(g))


def mul(f, g):
    """Multiplication of two polynomials via FFT convolution.

    Costs O(n log n) instead of O(n^2) naive multiplication.
    """
    return ifft(mul_fft(fft(f), fft(g)))


def div(f, g):
    """Division of two polynomials via FFT (coefficient representation)."""
    return ifft(div_fft(fft(f), fft(g)))


def adj(f):
    """Adjoint (conjugate in ring sense) of a polynomial (coefficient form)."""
    return ifft(adj_fft(fft(f)))


# ── Polynomial arithmetic (FFT representation) ───────────────────────────────

def add_fft(f_fft, g_fft):
    """Pointwise addition in FFT domain."""
    return add(f_fft, g_fft)


def sub_fft(f_fft, g_fft):
    """Pointwise subtraction in FFT domain."""
    return sub(f_fft, g_fft)


def mul_fft(f_fft, g_fft):
    """Pointwise multiplication in FFT domain (O(n) after transform)."""
    return [f_fft[i] * g_fft[i] for i in range(len(f_fft))]


def div_fft(f_fft, g_fft):
    """Pointwise division in FFT domain."""
    assert len(f_fft) == len(g_fft)
    return [f_fft[i] / g_fft[i] for i in range(len(f_fft))]


def adj_fft(f_fft):
    """Adjoint in FFT domain: pointwise complex conjugate."""
    return [c.conjugate() for c in f_fft]


# Ratio between ring degree and number of stored FFT coefficients (always 1 here)
fft_ratio = 1
