"""
Deterministic pseudo-random byte generator for the QBastion signer.

Wraps a ChaCha20 stream cipher seeded from a 56-byte signing seed. This
is used by the Falcon signing procedure to produce deterministic, reproducible
signatures when a seed is provided (as required by the NIST KAT vectors).

The generator exactly matches the interleaving strategy of the Falcon reference
C implementation: it produces 8 concurrent block-function calls at a time and
interleaves their 32-bit output words to match the reference byte ordering.

ChaCha20 matrix layout (column-major, 4x4 words):
    MAGIC[0]  MAGIC[1]  MAGIC[2]  MAGIC[3]
    seed[0]   seed[1]   seed[2]   seed[3]
    seed[4]   seed[5]   seed[6]   seed[7]
    seed[8]   seed[9]   ctr_lo    ctr_hi

Counter is a 64-bit integer split across words 12 and 13 of the seed.
"""

# Standard ChaCha20 magic constant words ("expand 32-byte k")
_CHACHA_MAGIC = [0x61707865, 0x3320646e, 0x79622d32, 0x6b206574]

# Keep legacy alias
CW = _CHACHA_MAGIC


def _rotate_left_32(word, shift):
    """Left-rotate a 32-bit integer by `shift` positions."""
    return ((word << shift) & 0xFFFFFFFF) | (word >> (32 - shift))


roll = _rotate_left_32  # legacy alias


class ChaCha20:
    """ChaCha20-based pseudo-random generator seeded from a 56-byte seed.

    Provides a `randombytes(k)` method that returns k pseudo-random bytes,
    matching the byte order used by the Falcon reference implementation's
    PRNG. The generator maintains an internal ring buffer of pre-generated
    hex bytes to amortize the cost of the block function.
    """

    def __init__(self, seed_bytes):
        """Initialize the generator from a 56-byte seed.

        The seed is parsed as 14 little-endian 32-bit words. Words 12 and 13
        are used as the low and high halves of a 64-bit counter.

        Args:
            seed_bytes: bytes-like object of length 56.
        """
        self._words = [
            int.from_bytes(seed_bytes[4 * i: 4 * (i + 1)], "little")
            for i in range(14)
        ]
        self._counter = self._words[12] + (self._words[13] << 32)
        self._hex_buf = ""

    def __repr__(self):
        parts = ', '.join(f'0x{w:08x}' for w in self._words)
        return f"ChaCha20(words=[{parts}], ctr={self._counter})"

    def _quarter_round(self, idx_a, idx_b, idx_c, idx_d):
        """Apply one ChaCha20 quarter-round to self._state in-place."""
        a = self._state[idx_a]
        b = self._state[idx_b]
        c = self._state[idx_c]
        d = self._state[idx_d]

        a = (a + b) & 0xFFFFFFFF;  d = _rotate_left_32(d ^ a, 16)
        c = (c + d) & 0xFFFFFFFF;  b = _rotate_left_32(b ^ c, 12)
        a = (a + b) & 0xFFFFFFFF;  d = _rotate_left_32(d ^ a, 8)
        c = (c + d) & 0xFFFFFFFF;  b = _rotate_left_32(b ^ c, 7)

        self._state[idx_a] = a
        self._state[idx_b] = b
        self._state[idx_c] = c
        self._state[idx_d] = d

    # Legacy method name expected by the block_update / test code
    qround = _quarter_round

    def update(self):
        """Run one ChaCha20 block function and return the 16-word output."""
        self._state = [0] * 16
        self._state[0:4]  = _CHACHA_MAGIC[:]
        self._state[4:14] = self._words[:10]
        self._state[14]   = self._words[10] ^ (self._counter & 0xFFFFFFFF)
        self._state[15]   = self._words[11] ^ (self._counter >> 32)
        snap = self._state[:]
        # 20 rounds = 10 double rounds
        for _ in range(10):
            # Column rounds
            self._quarter_round(0, 4, 8,  12)
            self._quarter_round(1, 5, 9,  13)
            self._quarter_round(2, 6, 10, 14)
            self._quarter_round(3, 7, 11, 15)
            # Diagonal rounds
            self._quarter_round(0, 5, 10, 15)
            self._quarter_round(1, 6, 11, 12)
            self._quarter_round(2, 7, 8,  13)
            self._quarter_round(3, 4, 9,  14)
        for i in range(16):
            self._state[i] = (self._state[i] + snap[i]) & 0xFFFFFFFF
        self._counter += 1
        return self._state

    # Expose internal state as `state` for backward compatibility
    @property
    def state(self):
        return self._state

    @state.setter
    def state(self, value):
        self._state = value

    # Expose _words as `s` for backward compatibility
    @property
    def s(self):
        return self._words

    @s.setter
    def s(self, value):
        self._words = value

    # Expose _counter as `ctr` for backward compat
    @property
    def ctr(self):
        return self._counter

    @ctr.setter
    def ctr(self, value):
        self._counter = value

    # Expose _hex_buf as `hexbytes` for backward compat
    @property
    def hexbytes(self):
        return self._hex_buf

    @hexbytes.setter
    def hexbytes(self, value):
        self._hex_buf = value

    def block_update(self):
        """Generate 8 parallel ChaCha20 blocks and interleave 32-bit words.

        This matches the reference code's strategy of processing 8 blocks
        simultaneously with outputs interleaved to improve SIMD utilization.

        Returns:
            A hex string representing 4096 pseudo-random bits.
        """
        interleaved = [None] * (16 * 8)
        for i in range(8):
            interleaved[i::8] = self.update()
        return "".join(w.to_bytes(4, "little").hex() for w in interleaved)

    def randombytes(self, k):
        """Return k pseudo-random bytes, refilling the buffer as needed.

        The byte ordering matches the Falcon reference implementation's
        interleaving convention.

        Args:
            k: number of bytes to generate

        Returns:
            bytes of length k.
        """
        if 2 * k > len(self._hex_buf):
            self._hex_buf = self.block_update()
        chunk = self._hex_buf[: 2 * k]
        # Byte-reverse each hex pair to match the reference endianness
        chunk = "".join(chunk[i: i + 2] for i in range(2 * k - 2, -1, -2))
        self._hex_buf = self._hex_buf[2 * k:]
        return bytes.fromhex(chunk)[::-1]
