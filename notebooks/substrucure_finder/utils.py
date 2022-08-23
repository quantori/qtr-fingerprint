from __future__ import annotations

from typing import List
import numpy as np
from pathlib import Path

from substrucure_finder.fingerprint import BitFingerprint, ByteFingerprint
from substrucure_finder import consts


def russelrao_radius(fingerprint: BitFingerprint) -> float:
    assert fingerprint.dtype == bool
    return 1 - np.sum(fingerprint) / len(fingerprint) + consts.radius_eps


def take_columns_from_bit_fingerprint(columns: List[int], fingerprint: BitFingerprint) -> BitFingerprint:
    assert fingerprint.dtype == bool
    return BitFingerprint(np.array([fingerprint[i] for i in columns]))


def is_sub_fingerprint(sub_fingerprint: np.ndarray, meta_fingerprint: np.ndarray) -> bool:
    assert sub_fingerprint.dtype == meta_fingerprint.dtype
    assert len(sub_fingerprint) == len(meta_fingerprint)
    return np.all(sub_fingerprint & meta_fingerprint == sub_fingerprint)


def byte_to_bits(byte: int) -> str:
    assert 0 <= byte < 256
    return bin(byte)[2::].rjust(8, '0')[::-1]


def bits_to_byte(bits: np.ndarray) -> int:
    assert 0 < len(bits) <= 8
    val = 0
    for i in range(len(bits)):
        assert 0 <= int(bits[i]) <= 1
        val |= int(bits[i]) << i
    val <<= 8 - len(bits)
    return val


def load_columns_from_file(file_path: Path) -> List[int]:
    assert file_path.is_file(), f"Path to load columns from must be a file, {file_path}"
    with Path(file_path).open('r') as f:
        return list(map(int, f.read().split()))


def save_columns_to_file(columns: List[int], file_path: Path) -> None:
    with Path(file_path).open('w') as f:
        f.write(' '.join(map(str, columns)))


def bit_fingerprint_to_str(fingerprint: BitFingerprint) -> str:
    assert fingerprint.dtype == bool
    return ''.join(map(lambda x: str(int(x)), fingerprint))


def bit_fingerprint_to_bytes(fingerprint: BitFingerprint) -> ByteFingerprint:
    assert fingerprint.dtype == bool
    values = [0] * ((len(fingerprint) + 7) >> 3)
    for i in range(0, len(fingerprint), 8):
        bits = fingerprint[i:i + 8]
        values[i >> 3] = bits_to_byte(bits)
    return ByteFingerprint(np.fromiter(values, dtype=np.uint8))


def byte_fingerprint_to_bits(fingerprint: ByteFingerprint, bits_cnt: int) -> BitFingerprint:
    assert fingerprint.dtype == np.uint8
    assert len(fingerprint) * 8 - 8 < bits_cnt <= len(fingerprint) * 8
    bits_arr = list(map(byte_to_bits, fingerprint))
    bits_arr[-1] = bits_arr[-1][(8 - bits_cnt % 8) % 8:]
    fp_str = ''.join(bits_arr)
    assert len(fp_str) == bits_cnt
    bit_fp = BitFingerprint(np.fromiter(map(int, fp_str), dtype=bool))
    return bit_fp
