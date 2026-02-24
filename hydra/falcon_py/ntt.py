"""
Number Theoretic Transform (NTT) over Z_q[x] / (x^n + 1) for QBastion.

This module performs fast polynomial arithmetic in the Falcon ring Z_q[x]/(phi),
where q = 12289 and phi = x^n + 1 for n a power of two.

The structure mirrors fft.py but works modulo the integer prime q rather than
over the complex numbers. All arithmetic is performed in constant-size integers
to avoid floating-point precision issues.
"""
from common import poly_split, poly_merge, FALCON_MODULUS

# Legacy aliases used throughout the project
split = poly_split
merge = poly_merge
q     = FALCON_MODULUS

# Import NTT-specific root tables and modular inverses
from ntt_constants import roots_dict_Zq, inv_mod_q


# ── Precomputed constants ─────────────────────────────────────────────────────

# Inverse of 2 modulo q (used in the radix-2 inverse NTT butterfly)
_INV2_MOD_Q = 6145
i2 = _INV2_MOD_Q  # legacy alias

# A primitive square root of -1 modulo q (i.e., i such that i^2 ≡ -1 mod q)
_SQRT_NEG1 = roots_dict_Zq[2][0]
sqr1 = _SQRT_NEG1  # legacy alias


# ── NTT butterfly helpers ─────────────────────────────────────────────────────

def _ntt_split(f_ntt):
    """Split an NTT polynomial into two half-degree NTT sub-polynomials.

    Uses the NTT analogue of splitfft_2 over Z_q.

    Args:
        f_ntt: polynomial in NTT representation (length n)

    Returns:
        [f0_ntt, f1_ntt] each of length n/2, both in NTT form
    """
    deg   = len(f_ntt)
    roots = roots_dict_Zq[deg]
    half  = deg // 2
    f0 = [0] * half
    f1 = [0] * half
    for k in range(half):
        lo = f_ntt[2 * k]
        hi = f_ntt[2 * k + 1]
        f0[k] = (_INV2_MOD_Q * (lo + hi)) % q
        f1[k] = (_INV2_MOD_Q * (lo - hi) * inv_mod_q[roots[2 * k]]) % q
    return [f0, f1]


split_ntt = _ntt_split  # public alias


def _ntt_merge(pair_ntt):
    """Merge two half-degree NTT polynomials into a full-degree NTT polynomial.

    Args:
        pair_ntt: [f0_ntt, f1_ntt] each of length n/2

    Returns:
        The merged NTT polynomial of length n.
    """
    f0_ntt, f1_ntt = pair_ntt
    half   = len(f0_ntt)
    deg    = 2 * half
    roots  = roots_dict_Zq[deg]
    result = [0] * deg
    for k in range(half):
        w_f1 = roots[2 * k] * f1_ntt[k]
        result[2 * k]     = (f0_ntt[k] + w_f1) % q
        result[2 * k + 1] = (f0_ntt[k] - w_f1) % q
    return result


merge_ntt = _ntt_merge  # public alias


# ── Forward and inverse NTT ───────────────────────────────────────────────────

def ntt(poly):
    """Compute the NTT of a polynomial over Z_q.

    Recursive split-radix transform, analogous to fft() but over the finite
    field Z_q rather than the complex numbers.

    Args:
        poly: coefficient list of length n (values in Z_q)

    Returns:
        NTT representation of length n.
    """
    deg = len(poly)
    if deg > 2:
        even, odd     = poly_split(poly)
        even_ntt      = ntt(even)
        odd_ntt       = ntt(odd)
        return _ntt_merge([even_ntt, odd_ntt])
    # Base case: 2-point DFT over Z_q using the primitive sqrt(-1)
    return [
        (poly[0] + _SQRT_NEG1 * poly[1]) % q,
        (poly[0] - _SQRT_NEG1 * poly[1]) % q,
    ]


def intt(f_ntt):
    """Compute the inverse NTT, recovering a polynomial over Z_q.

    Args:
        f_ntt: NTT representation (list of length n, values in Z_q)

    Returns:
        Coefficient list of length n (values in Z_q).
    """
    deg = len(f_ntt)
    if deg > 2:
        f0_ntt, f1_ntt = _ntt_split(f_ntt)
        even           = intt(f0_ntt)
        odd            = intt(f1_ntt)
        return poly_merge([even, odd])
    # Base case
    return [
        (_INV2_MOD_Q * (f_ntt[0] + f_ntt[1])) % q,
        (_INV2_MOD_Q * inv_mod_q[_SQRT_NEG1] * (f_ntt[0] - f_ntt[1])) % q,
    ]


# ── Polynomial arithmetic over Z_q (coefficient form) ────────────────────────

def add_zq(f, g):
    """Coefficient-wise addition modulo q."""
    assert len(f) == len(g)
    return [(f[i] + g[i]) % q for i in range(len(f))]


def neg_zq(f):
    """Coefficient-wise negation modulo q."""
    return [(-c) % q for c in f]


def sub_zq(f, g):
    """Coefficient-wise subtraction modulo q."""
    return add_zq(f, neg_zq(g))


def mul_zq(f, g):
    """Polynomial multiplication modulo (x^n + 1, q) via NTT convolution."""
    return intt(mul_ntt(ntt(f), ntt(g)))


def div_zq(f, g):
    """Polynomial division modulo (x^n + 1, q) via NTT.

    Raises ZeroDivisionError if g is not invertible in the NTT domain.
    """
    return intt(div_ntt(ntt(f), ntt(g)))


# ── Polynomial arithmetic (NTT representation) ───────────────────────────────

def add_ntt(f_ntt, g_ntt):
    """Pointwise addition in the NTT domain (same as add_zq)."""
    return add_zq(f_ntt, g_ntt)


def sub_ntt(f_ntt, g_ntt):
    """Pointwise subtraction in the NTT domain."""
    return sub_zq(f_ntt, g_ntt)


def mul_ntt(f_ntt, g_ntt):
    """Pointwise multiplication in the NTT domain, modulo q."""
    assert len(f_ntt) == len(g_ntt)
    return [(f_ntt[i] * g_ntt[i]) % q for i in range(len(f_ntt))]


def div_ntt(f_ntt, g_ntt):
    """Pointwise division in the NTT domain using precomputed inverses.

    Raises ZeroDivisionError if any entry of g_ntt is zero modulo q.
    """
    assert len(f_ntt) == len(g_ntt)
    if any(v == 0 for v in g_ntt):
        raise ZeroDivisionError("Cannot invert a zero NTT coefficient")
    return [(f_ntt[i] * inv_mod_q[g_ntt[i]]) % q for i in range(len(f_ntt))]


# Ratio between ring degree n and the number of stored NTT coefficients
ntt_ratio = 1
