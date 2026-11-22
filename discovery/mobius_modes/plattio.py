"""Low-level reader for Dave Platt's zero tables (format read off Jonathan Bober's platt_zeros.py).
File: 8-byte block count; each block: 32-byte header (t0, t1: double; Nt0, Nt1: uint64) then (Nt1-Nt0) entries of
13 bytes (uint64 z1, uint32 z2, uint8 z3).  Each entry is the integer gap, in units 2**-101, from the previous
zero (the first entry of a block: from t0).  Zero index N is 1-based: N = Nt0 + i + 1 for entry i."""
import struct
import numpy as np
from fractions import Fraction

DT13 = np.dtype([('a', '<u8'), ('b', '<u4'), ('c', 'u1')])
assert DT13.itemsize == 13
EPS_EXP = 101


def blocks(path):
    """Yield (block_number, offset, t0, t1, Nt0, Nt1) for every block of a file, plus the file size check."""
    with open(path, 'rb') as f:
        nb = struct.unpack('Q', f.read(8))[0]
        off = 8
        for b in range(nb):
            f.seek(off)
            hdr = f.read(32)
            if len(hdr) < 32:
                raise ValueError(f"{path}: truncated header at block {b}")
            t0, t1, N0, N1 = struct.unpack('ddQQ', hdr)
            if N1 < N0:
                raise ValueError(f"{path}: bad block {b}")
            yield b, off, t0, t1, N0, N1
            off += 32 + 13 * (N1 - N0)


def read_block(path, offset, n):
    with open(path, 'rb') as f:
        f.seek(offset + 32)
        return np.frombuffer(f.read(13 * n), dtype=DT13)


def delta_float(arr):
    """Entries as float64 gaps (exact to 1e-16 relative: built from the integer pieces)."""
    return (arr['a'].astype(np.float64) * 2.0 ** -101 + arr['b'].astype(np.float64) * 2.0 ** -37
            + arr['c'].astype(np.float64) * 2.0 ** -5)


def delta_sum_int(arr):
    """Exact Python-int sum of the integer entries, in units 2**-101."""
    a = arr['a']
    lo = int((a & np.uint64(0xFFFFFFFF)).sum(dtype=np.uint64))
    hi = int((a >> np.uint64(32)).sum(dtype=np.uint64))
    b = int(arr['b'].sum(dtype=np.uint64))
    c = int(arr['c'].sum(dtype=np.uint64))
    return lo + (hi << 32) + (b << 64) + (c << 96)


def t0_int(t0):
    """t0 (a double) as an exact integer in units 2**-101."""
    fr = Fraction(t0) * (1 << EPS_EXP)
    assert fr.denominator == 1
    return int(fr)


def entry_int(arr, i):
    return int(arr['a'][i]) + (int(arr['b'][i]) << 64) + (int(arr['c'][i]) << 96)
