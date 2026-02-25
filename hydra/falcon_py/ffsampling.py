"""
Fast Fourier lattice algorithms for the QBastion signing backend.

Implements three key subroutines of the Falcon signing procedure:

1. **Gram** — Computes the Gram matrix G = B * B^* for a basis B.
2. **ffLDL** — Recursively decomposes the Gram matrix into an LDL^* tree
   (a binary tree of polynomials) using the fast Fourier lifting trick, as
   described in Algorithm 9 of the Falcon specification.
3. **ffSampling** — Samples a short lattice vector near a target using the
   ffNP (Fast Fourier nearest-plane) algorithm (Algorithm 11 of Falcon).

All functions are available in both coefficient and FFT representations.
FFT variants are preferred in production for their O(n log n) complexity.
"""
from common import split, merge
from fft import add, sub, mul, div, adj
from fft import add_fft, sub_fft, mul_fft, div_fft, adj_fft
from fft import split_fft, merge_fft, fft_ratio
from samplerz import samplerz


# ── Gram matrix ───────────────────────────────────────────────────────────────

def gram(basis):
    """Compute the Gram matrix G = B * B^* for a basis matrix B.

    For a 2x2 matrix B with polynomial entries, computes the Hermitian
    positive-definite Gram matrix G[i][j] = sum_k B[i][k] * adj(B[j][k]).

    Args:
        basis: a 2x2 list of polynomial lists, in coefficient representation.

    Returns:
        G: a 2x2 Gram matrix of polynomials (coefficient representation).
    """
    row_indices = range(len(basis))
    num_cols    = len(basis[0])
    degree      = len(basis[0][0])
    G = [[[0] * degree for _ in row_indices] for _ in row_indices]
    for i in row_indices:
        for j in row_indices:
            for k in range(num_cols):
                G[i][j] = add(G[i][j], mul(basis[i][k], adj(basis[j][k])))
    return G


# ── LDL decomposition ─────────────────────────────────────────────────────────

def ldl(G):
    """Compute the LDL* decomposition of a 2x2 polynomial Gram matrix G.

    Produces lower triangular L and diagonal D such that G = L * D * L^*.
    This is the polynomial analogue of Algorithm 8 from the Falcon spec.

    Args:
        G: a 2x2 Gram matrix in coefficient representation.

    Returns:
        [L, D]: both are 2x2 polynomial matrices (coefficient representation).
    """
    degree = len(G[0][0])
    assert len(G) == 2 and len(G[0]) == 2
    zero = [0] * degree
    one  = [1] + [0] * (degree - 1)
    d00  = G[0][0][:]
    l10  = div(G[1][0], G[0][0])
    d11  = sub(G[1][1], mul(mul(l10, adj(l10)), G[0][0]))
    L = [[one, zero], [l10, one]]
    D = [[d00, zero], [zero, d11]]
    return [L, D]


def ldl_fft(G):
    """LDL* decomposition of a 2x2 Gram matrix in FFT representation.

    Equivalent to ldl() but operates on polynomials in FFT form, avoiding
    the need to convert back and forth for each operation.

    Args:
        G: a 2x2 Gram matrix in FFT representation.

    Returns:
        [L, D]: both are 2x2 polynomial matrices (FFT representation).
    """
    degree = len(G[0][0])
    assert len(G) == 2 and len(G[0]) == 2
    zero = [0] * degree
    one  = [1] * degree
    d00  = G[0][0][:]
    l10  = div_fft(G[1][0], G[0][0])
    d11  = sub_fft(G[1][1], mul_fft(mul_fft(l10, adj_fft(l10)), G[0][0]))
    L = [[one, zero], [l10, one]]
    D = [[d00, zero], [zero, d11]]
    return [L, D]


# ── Fast Fourier LDL tree (ffLDL) ────────────────────────────────────────────

def ffldl(G):
    """Recursively build the ffLDL* decomposition tree of a Gram matrix.

    Implements Algorithm 9 (ffLDL) from the Falcon specification in the
    coefficient representation. The tree structure enables O(n log n)
    ffSampling — each internal node stores an L10 entry, and leaves store
    the diagonal scalars.

    Args:
        G: a 2x2 Gram matrix (coefficient representation).

    Returns:
        A tree [L10, subtree_0, subtree_1] where subtrees recurse to leaves
        [L10, d00, d11] (coefficient representation, degree 2).
    """
    deg = len(G[0][0])
    L, D = ldl(G)
    if deg > 2:
        d00, d01 = split(D[0][0])
        d10, d11 = split(D[1][1])
        G0 = [[d00, d01], [adj(d01), d00]]
        G1 = [[d10, d11], [adj(d11), d10]]
        return [L[1][0], ffldl(G0), ffldl(G1)]
    # Base case: degree-2 polynomials — only the constant term matters
    D[0][0][1] = 0
    D[1][1][1] = 0
    return [L[1][0], D[0][0], D[1][1]]


def ffldl_fft(G):
    """Build the ffLDL* tree in FFT representation (production path).

    Same recursive structure as ffldl() but uses FFT-domain arithmetic
    throughout, giving O(n log n) overall complexity for tree construction.

    Args:
        G: a 2x2 Gram matrix in FFT representation.

    Returns:
        A tree [L10_fft, subtree_0, subtree_1], leaves are single-element
        FFT lists [sigma^2] (one complex value each, purely real).
    """
    deg = len(G[0][0]) * fft_ratio
    L, D = ldl_fft(G)
    if deg > 2:
        d00, d01 = split_fft(D[0][0])
        d10, d11 = split_fft(D[1][1])
        G0 = [[d00, d01], [adj_fft(d01), d00]]
        G1 = [[d10, d11], [adj_fft(d11), d10]]
        return [L[1][0], ffldl_fft(G0), ffldl_fft(G1)]
    # Base case: each element is a single complex value (real part = variance)
    return [L[1][0], D[0][0], D[1][1]]


# ── Fast Fourier nearest-plane (ffNP) ─────────────────────────────────────────

def ffnp(target, ldl_tree):
    """Run the ffNP nearest-plane reducer in coefficient representation.

    Given a target vector and an ffLDL* tree, produce the lattice point
    closest to target using the recursive fast-Fourier nearest-plane heuristic.
    Used primarily for testing; the FFT version is faster in practice.

    Args:
        target:   a 2-vector of polynomials [t0, t1] (coefficient form).
        ldl_tree: the ffLDL* tree for this lattice basis.

    Returns:
        z: a 2-vector of integer-rounded polynomials forming the reduced point.
    """
    deg = len(target[0])
    z   = [None, None]
    if deg > 1:
        l10, T0, T1 = ldl_tree
        z[1]  = merge(ffnp(split(target[1]), T1))
        t0b   = add(target[0], mul(sub(target[1], z[1]), l10))
        z[0]  = merge(ffnp(split(t0b), T0))
    else:
        # Base case: round each scalar to the nearest integer
        z[0] = [round(target[0][0])]
        z[1] = [round(target[1][0])]
    return z


def ffnp_fft(target_fft, ldl_tree):
    """Run the ffNP nearest-plane reducer in FFT representation.

    The FFT variant of ffnp(): avoids intermediate conversions by working
    entirely in the FFT domain. The base case reads the real part of each
    FFT coefficient (which equals the single polynomial coefficient at n=1).

    Args:
        target_fft: a 2-vector [t0_fft, t1_fft] (FFT representation).
        ldl_tree:   the ffLDL* tree built by ffldl_fft().

    Returns:
        z: a 2-vector of integer-rounded polynomials (FFT form, length 1 at base).
    """
    deg = len(target_fft[0]) * fft_ratio
    z   = [0, 0]
    if deg > 1:
        l10, T0, T1 = ldl_tree
        z[1]  = merge_fft(ffnp_fft(split_fft(target_fft[1]), T1))
        t0b   = add_fft(target_fft[0], mul_fft(sub_fft(target_fft[1], z[1]), l10))
        z[0]  = merge_fft(ffnp_fft(split_fft(t0b), T0))
    else:
        z[0] = [round(target_fft[0][0].real)]
        z[1] = [round(target_fft[1][0].real)]
    return z


# ── Fast Fourier sampling (ffSampling) ────────────────────────────────────────

def ffsampling_fft(target_fft, ldl_tree, sigmin, randombytes):
    """Sample a lattice vector close to `target_fft` using Gaussian sampling.

    Implements Algorithm 11 (ffSampling) from the Falcon specification.
    Replaces the deterministic rounding of ffNP with a discrete Gaussian
    sampler at each leaf, producing a randomized short preimage vector.

    Args:
        target_fft:  a 2-vector [t0_fft, t1_fft] in FFT representation.
        ldl_tree:    the normalized ffLDL* tree (leaves store sigma / ||b_i||).
        sigmin:      per-coordinate Gaussian floor (from FalconParam.sigmin).
        randombytes: entropy source callable, signature `randombytes(k) -> bytes`.

    Returns:
        z: a 2-vector of sampled integer polynomials (FFT representation).
    """
    deg = len(target_fft[0]) * fft_ratio
    z   = [0, 0]
    if deg > 1:
        l10, T0, T1 = ldl_tree
        z[1]  = merge_fft(ffsampling_fft(split_fft(target_fft[1]), T1, sigmin, randombytes))
        t0b   = add_fft(target_fft[0], mul_fft(sub_fft(target_fft[1], z[1]), l10))
        z[0]  = merge_fft(ffsampling_fft(split_fft(t0b), T0, sigmin, randombytes))
    else:
        # Leaf node: sample from D_{Z, t_i.real, T[0]}
        leaf_sigma = ldl_tree[0]
        z[0] = [samplerz(target_fft[0][0].real, leaf_sigma, sigmin, randombytes)]
        z[1] = [samplerz(target_fft[1][0].real, leaf_sigma, sigmin, randombytes)]
    return z
