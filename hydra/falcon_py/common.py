"""
Shared utilities and constants for the QBastion Falcon integration layer.
This module provides core polynomial manipulation routines used throughout
the post-quantum signature pipeline.
"""


# Falcon modulus: q = 12 * 2^10 + 1
FALCON_MODULUS = 12 * 1024 + 1
# Alias for interoperability with other modules that import `q` directly
q = FALCON_MODULUS


def poly_split(poly):
    """Decompose a polynomial into two interleaved sub-polynomials.

    Given poly = [a0, a1, a2, a3, ...], produces:
        even = [a0, a2, ...]
        odd  = [a1, a3, ...]

    Args:
        poly: input polynomial in coefficient representation

    Returns:
        A pair [even_coeffs, odd_coeffs]
    """
    length = len(poly)
    half = length // 2
    even_coeffs = [poly[2 * idx] for idx in range(half)]
    odd_coeffs  = [poly[2 * idx + 1] for idx in range(half)]
    return [even_coeffs, odd_coeffs]


# Legacy alias used by fft.py / ntt.py
split = poly_split


def poly_merge(sub_polys):
    """Reconstruct a polynomial from two interleaved sub-polynomials.

    The inverse of poly_split: interleaves even and odd coefficient lists
    back into a single polynomial.

    Args:
        sub_polys: a list [even_coeffs, odd_coeffs]

    Returns:
        The reconstructed polynomial
    """
    even_part, odd_part = sub_polys
    total_len = 2 * len(even_part)
    result = [0] * total_len
    for idx in range(total_len // 2):
        result[2 * idx]     = even_part[idx]
        result[2 * idx + 1] = odd_part[idx]
    return result


# Legacy alias
merge = poly_merge


def squared_norm(vec):
    """Compute the squared Euclidean norm of a 2D polynomial vector.

    Args:
        vec: a list of polynomials, e.g. [s0, s1]

    Returns:
        The sum of squares of all coefficients across all component polynomials.
    """
    total = 0
    for component in vec:
        for coeff in component:
            total += coeff * coeff
    return total


# Legacy alias used in test.py
sqnorm = squared_norm
