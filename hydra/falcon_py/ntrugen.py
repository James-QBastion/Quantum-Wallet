"""
NTRU key generation for the QBastion Falcon signer.

Implements Section 3.8.2 of the Falcon specification:
- Karatsuba polynomial multiplication
- Field norm and Galois conjugate
- NTRU equation solver (NTRUSolve)
- Babai reduction (Reduce)
- Key pair generation (NTRUGen, Algorithm 5)

All polynomial types are list[int] or list[float] throughout.
The ring is Z[x] / (x^n + 1) for n a power of two.
"""
from __future__ import annotations

from fft import fft, ifft, add_fft, mul_fft, adj_fft, div_fft
from fft import add, mul, div, adj
from ntt import ntt
from common import squared_norm as sqnorm
from samplerz import samplerz

# Alias sqnorm for use in gs_norm
_sqnorm = sqnorm

# Falcon modulus (same as common.FALCON_MODULUS)
q: int = 12 * 1024 + 1


# ── Polynomial multiplication ─────────────────────────────────────────────────

def karatsuba(a: list, b: list, n: int) -> list:
    """Schoolbook Karatsuba multiplication of two degree-n polynomials.

    Produces a degree-2n polynomial ab = a * b whose coefficients may be
    integers or floats. Reduction modulo (x^n + 1) is NOT performed here —
    see karamul() for the reduced version.

    Args:
        a: coefficient list of length n
        b: coefficient list of length n
        n: degree (must be a power of two)

    Returns:
        ab: coefficient list of length 2n
    """
    if n == 1:
        return [a[0] * b[0], 0]

    half = n // 2
    a0, a1 = a[:half], a[half:]
    b0, b1 = b[:half], b[half:]
    # Middle sum trick to save one recursive call
    sum_a = [a0[i] + a1[i] for i in range(half)]
    sum_b = [b0[i] + b1[i] for i in range(half)]

    lo  = karatsuba(a0, b0, half)
    hi  = karatsuba(a1, b1, half)
    mid = karatsuba(sum_a, sum_b, half)

    for i in range(n):
        mid[i] -= lo[i] + hi[i]

    ab = [0] * (2 * n)
    for i in range(n):
        ab[i]        += lo[i]
        ab[i + n]    += hi[i]
        ab[i + half] += mid[i]
    return ab


def karamul(a: list[int], b: list[int]) -> list[int]:
    """Karatsuba multiplication reduced modulo (x^n + 1).

    Args:
        a: polynomial of degree < n
        b: polynomial of degree < n

    Returns:
        a * b mod (x^n + 1), a list of n integers
    """
    n  = len(a)
    ab = karatsuba(a, b, n)
    return [ab[i] - ab[i + n] for i in range(n)]


# ── Ring-theoretic helpers ────────────────────────────────────────────────────

def galois_conjugate(a: list[int]) -> list[int]:
    """Compute the Galois conjugate of a ∈ Q[x] / (x^n + 1).

    The Galois conjugate of a(x) is a(-x), obtained by flipping the sign of
    every odd-indexed coefficient.

    Args:
        a: coefficient list of length n

    Returns:
        Galois conjugate as a coefficient list of length n.
    """
    return [((-1) ** i) * a[i] for i in range(len(a))]


def field_norm(a: list[int]) -> list[int]:
    """Project a ∈ Q[x]/(x^n+1) to Q[x]/(x^(n/2)+1) via the field norm.

    Computes N(a) = a_e^2 - x * a_o^2 where a_e (resp. a_o) are the
    coefficients of even (resp. odd) degree in a.

    Args:
        a: coefficient list of length n (n must be a power of two, n >= 2)

    Returns:
        Field norm N(a) as a coefficient list of length n/2.
    """
    half = len(a) // 2
    a_even = [a[2 * i]     for i in range(half)]
    a_odd  = [a[2 * i + 1] for i in range(half)]
    sq_even = karamul(a_even, a_even)
    sq_odd  = karamul(a_odd,  a_odd)
    result  = sq_even[:]
    for i in range(half - 1):
        result[i + 1] -= sq_odd[i]
    result[0] += sq_odd[half - 1]
    return result


def lift(a: list[int]) -> list[int]:
    """Lift a ∈ Q[x]/(x^(n/2)+1) to Q[x]/(x^n+1) by substituting x → x^2.

    Args:
        a: coefficient list of length n/2

    Returns:
        a(x^2) as a coefficient list of length n (odd-indexed entries are 0).
    """
    n      = len(a)
    result = [0] * (2 * n)
    for i in range(n):
        result[2 * i] = a[i]
    return result


# ── Arithmetic helpers ────────────────────────────────────────────────────────

def bitsize(a: int) -> int:
    """Compute the byte-rounded bit-length of |a| (ignoring sign).

    Returns the bit-length rounded up to the next multiple of 8. Used to
    adaptively scale polynomial coefficients in the Babai reduction.

    Args:
        a: an integer (may be negative)

    Returns:
        Bit size rounded to next multiple of 8.
    """
    val, bits = abs(a), 0
    while val:
        bits += 8
        val >>= 8
    return bits


# ── Babai reduction ───────────────────────────────────────────────────────────

def reduce(
    f: list[int],
    g: list[int],
    F: list[int],
    G: list[int],
) -> tuple[list[int], list[int]]:
    """Perform Babai's nearest-plane reduction of (F, G) with respect to (f, g).

    Updates (F, G) ← (F, G) − k*(f, g) where k = round((F*f̄ + G*ḡ) / (f*f̄ + g*ḡ)).
    This is Algorithm 7 (Reduce) from the Falcon specification.

    Multiple iterations may be needed when working in finite precision over
    very large polynomial coefficients.

    Args:
        f, g: short NTRU basis polynomials
        F, G: tall NTRU output polynomials (modified in place)

    Returns:
        (F, G) after reduction.
    """
    n    = len(f)
    size = max(53,
               bitsize(min(f)), bitsize(max(f)),
               bitsize(min(g)), bitsize(max(g)))

    f_sc   = [c >> (size - 53) for c in f]
    g_sc   = [c >> (size - 53) for c in g]
    fa_fft = fft(f_sc)
    ga_fft = fft(g_sc)

    while True:
        big = max(53,
                  bitsize(min(F)), bitsize(max(F)),
                  bitsize(min(G)), bitsize(max(G)))
        if big < size:
            break
        F_sc   = [c >> (big - 53) for c in F]
        G_sc   = [c >> (big - 53) for c in G]
        Fa_fft = fft(F_sc)
        Ga_fft = fft(G_sc)

        den_fft = add_fft(mul_fft(fa_fft, adj_fft(fa_fft)),
                          mul_fft(ga_fft, adj_fft(ga_fft)))
        num_fft = add_fft(mul_fft(Fa_fft, adj_fft(fa_fft)),
                          mul_fft(Ga_fft, adj_fft(ga_fft)))
        k_fft = div_fft(num_fft, den_fft)
        k     = [int(round(c)) for c in ifft(k_fft)]
        if all(c == 0 for c in k):
            break
        fk = karamul(f, k)
        gk = karamul(g, k)
        shift = big - size
        for i in range(n):
            F[i] -= fk[i] << shift
            G[i] -= gk[i] << shift
    return F, G


# ── Extended GCD ──────────────────────────────────────────────────────────────

def xgcd(b: int, n: int) -> tuple[int, int, int]:
    """Extended Euclidean algorithm for integers b and n.

    Returns (d, u, v) such that d = gcd(b, n) = u*b + v*n.

    Args:
        b, n: non-negative integers

    Returns:
        (d, u, v): gcd and Bézout coefficients.
    """
    x0, x1, y0, y1 = 1, 0, 0, 1
    while n:
        q, b, n = b // n, n, b % n
        x0, x1  = x1, x0 - q * x1
        y0, y1  = y1, y0 - q * y1
    return b, x0, y0


# ── NTRU equation solver ──────────────────────────────────────────────────────

def ntru_solve(f: list[int], g: list[int]) -> tuple[list[int], list[int]]:
    """Solve the NTRU equation: find F, G in Z[x]/(x^n+1) such that f*G - g*F = q.

    Implements NTRUSolve from the Falcon specification using a recursive
    field-norm descent. Returns (F, G) after Babai reduction.

    Args:
        f, g: NTRU key polynomials (short, invertible mod q)

    Returns:
        (F, G): solution polynomials

    Raises:
        ValueError: if gcd(f[0], g[0]) ≠ 1 at the base case.
    """
    n = len(f)
    if n == 1:
        d, u, v = xgcd(f[0], g[0])
        if d != 1:
            raise ValueError(f"Non-invertible base case: gcd({f[0]}, {g[0]}) = {d}")
        return [-q * v], [q * u]

    fp  = field_norm(f)
    gp  = field_norm(g)
    Fp, Gp = ntru_solve(fp, gp)
    F = karamul(lift(Fp), galois_conjugate(g))
    G = karamul(lift(Gp), galois_conjugate(f))
    return reduce(f, g, F, G)


# ── Gram-Schmidt norm ─────────────────────────────────────────────────────────

def gs_norm(f: list[int], g: list[int], modulus: int) -> float:
    """Compute the squared Gram-Schmidt norm of the NTRU lattice basis [[g, -f], [G, -F]].

    This is equivalent to line 9 of Algorithm 5 (NTRUGen) from the Falcon spec.
    Used to enforce the key quality bound σ_GS^2 ≤ 1.17^2 * q.

    Args:
        f, g:    NTRU short polynomials
        modulus: the Falcon prime q

    Returns:
        max(||[f, g]||^2, q^2 * ||[F*, G*]||^2) as a float.
    """
    sq_fg = _sqnorm([f, g])
    ffgg  = add(mul(f, adj(f)), mul(g, adj(g)))
    Ft    = div(adj(g), ffgg)
    Gt    = div(adj(f), ffgg)
    sq_FG = (modulus ** 2) * _sqnorm([Ft, Gt])
    return max(sq_fg, sq_FG)


# ── Polynomial generator ──────────────────────────────────────────────────────

def gen_poly(n: int) -> list[int]:
    """Sample a short polynomial of degree < n from D_{Z, 0, σ_fg}.

    Uses the additive Gaussian trick: sums 4096/n independent Gaussian samples
    to produce a sample from a Gaussian of std dev sqrt(4096/n) * σ_base, which
    matches σ_fg = 1.17 * sqrt(q / (2n)).

    Args:
        n: ring degree (power of two, < 4096)

    Returns:
        A list of n integers sampled from the target distribution.
    """
    # σ_fg ≈ 1.43300980528773 for all n (derived from 1.17 * sqrt(12289/8192))
    sigma_base = 1.43300980528773
    assert n < 4096

    oversample = 4096
    raw = [samplerz(0.0, sigma_base, sigma_base - 0.001) for _ in range(oversample)]
    chunk = oversample // n
    return [sum(raw[i * chunk + j] for j in range(chunk)) for i in range(n)]


# ── Main key generation ───────────────────────────────────────────────────────

def ntru_gen(n: int) -> tuple[list[int], list[int], list[int], list[int]]:
    """Generate an NTRU key pair (f, g, F, G) for Falcon-n.

    Implements Algorithm 5 (NTRUGen) from the Falcon specification. Repeatedly
    samples (f, g) until:
      - The Gram-Schmidt norm is within the 1.17^2 * q bound
      - f is invertible in Z_q (no zero NTT coefficient)
      - The NTRU equation f*G - g*F = q has an integer solution

    Args:
        n: ring degree (must be a power-of-two in [2, 1024])

    Returns:
        (f, g, F, G) — four polynomials in Z[x]/(x^n+1)
    """
    gs_bound = (1.17 ** 2) * q
    while True:
        f = gen_poly(n)
        g = gen_poly(n)
        if gs_norm(f, g, q) > gs_bound:
            continue
        if any(v == 0 for v in ntt(f)):
            continue
        try:
            F, G = ntru_solve(f, g)
            return f, g, [int(c) for c in F], [int(c) for c in G]
        except ValueError:
            continue
