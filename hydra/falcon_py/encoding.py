"""
Signature codec for the QBastion lattice signing pipeline.

Encodes and decodes polynomial signatures using a hybrid sign+unary bit packing
scheme, optimized for compact lattice signatures. This encoding is compatible
with the Falcon NIST submission specification.

Encoding layout per coefficient:
  [1 bit sign] [7 bits low] [unary-encoded high bits] [terminator '1'] 
"""


def pack_signature(coeffs, byte_budget):
    """Encode a list of integer coefficients into a fixed-length byte buffer.

    Uses a hybrid sign + lower-7-bits + unary-upper-bits scheme. If the
    total bitstring exceeds 8 * byte_budget bits, returns False to
    signal the retry loop to resample.

    Args:
        coeffs:      list of integers (the signature polynomial s1)
        byte_budget: target byte length of the output

    Returns:
        A bytes object of exactly `byte_budget` bytes, or False on overflow.
    """
    bit_buf = ""
    for val in coeffs:
        # Sign bit: '1' for negative, '0' for non-negative
        sign_bit = "1" if val < 0 else "0"
        abs_val  = abs(val)
        # Lower 7 bits in binary (fixed width)
        low_bits = format(abs_val % (1 << 7), '#09b')[2:]
        # High bits encoded in unary: k zeros then a terminating '1'
        unary_high = "0" * (abs_val >> 7) + "1"
        bit_buf += sign_bit + low_bits + unary_high

    if len(bit_buf) > 8 * byte_budget:
        return False

    # Zero-pad to fill the buffer
    bit_buf += "0" * (8 * byte_budget - len(bit_buf))
    octets = [int(bit_buf[8 * i: 8 * i + 8], 2) for i in range(len(bit_buf) // 8)]
    return bytes(octets)


# Public alias matching the original API contract
compress = pack_signature


def unpack_signature(raw_bytes, byte_budget, degree):
    """Decode a packed byte buffer back into a polynomial coefficient list.

    Args:
        raw_bytes:   the byte-encoded signature buffer
        byte_budget: expected byte length (checked for consistency)
        degree:      number of coefficients to recover (ring degree n)

    Returns:
        A list of `degree` integers, or False if the encoding is invalid.
    """
    if len(raw_bytes) > byte_budget:
        # Encoding too long — reject without printing to stderr
        return False

    # Convert each byte into 8 bits (without the leading `1` from bin())
    bit_str = ""
    for byte_val in raw_bytes:
        bit_str += bin((1 << 8) ^ byte_val)[3:]

    recovered = []

    # Strip trailing zero padding bits
    while bit_str and bit_str[-1] == "0":
        bit_str = bit_str[:-1]

    try:
        while bit_str and len(recovered) < degree:
            # Recover sign
            polarity = -1 if bit_str[0] == "1" else 1
            # Recover low 7 bits
            low = int(bit_str[1:8], 2)
            # Recover unary-encoded high bits
            cursor, high_part = 8, 0
            while bit_str[cursor] == "0":
                cursor += 1
                high_part += 1
            coeff = polarity * (low + (high_part << 7))
            # Enforce canonical zero encoding: -0 is invalid
            if coeff == 0 and polarity == -1:
                return False
            recovered.append(coeff)
            bit_str = bit_str[cursor + 1:]

        if len(recovered) != degree:
            return False
        return recovered
    except IndexError:
        return False


# Public alias
decompress = unpack_signature
