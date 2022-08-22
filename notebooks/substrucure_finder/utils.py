from __future__ import annotations

from typing import List
import numpy as np
from pathlib import Path

from substrucure_finder.fingerprint import BitFingerprint, ByteFingerprint
from substrucure_finder import consts


def russelrao_radius(fingerprint: BitFingerprint) -> float:
    return 1 - sum(fingerprint) / len(fingerprint) + consts.radius_eps


def take_columns_from_bit_fingerprint(columns: List[int], fingerprint: BitFingerprint) -> BitFingerprint:
    return BitFingerprint(np.array([fingerprint[i] for i in columns]))


def is_sub_fingerprint(sub_fingerprint: np.ndarray, meta_fingerprint: np.ndarray) -> bool:
    assert len(sub_fingerprint) == len(meta_fingerprint)
    return np.all(sub_fingerprint & meta_fingerprint == sub_fingerprint)


def byte_to_bits(byte: int) -> str:
    assert 0 <= byte < 256
    return bin(byte)[2::].rjust(8, '0')[::-1]


def load_columns_from_file(file_path: Path) -> List[int]:
    assert file_path.is_file(), f"Path to load columns from must be a file, {file_path}"
    with Path(file_path).open('r') as f:
        return list(map(int, f.read().split()))


def save_columns_to_file(columns: List[int], file_path: Path) -> None:
    with Path(file_path).open('w') as f:
        f.write(' '.join(map(str, columns)))


def bit_fingerprint_to_str(fingerprint: BitFingerprint) -> str:
    return ''.join(map(str, fingerprint))


def bit_fingerprint_to_bytes(fingerprint: BitFingerprint) -> ByteFingerprint:
    values = [0] * ((len(fingerprint) + 7) >> 3)
    for i in range(0, len(fingerprint), 8):
        val = 0
        for j in range(i, min(i + 8, len(fingerprint))):
            val = val * 2 + int(fingerprint[j])
        values[i >> 3] = val
    return ByteFingerprint(np.fromiter(values, dtype=np.uint8))


def byte_fingerprint_to_bits(fingerprint: ByteFingerprint, bits_cnt: int) -> BitFingerprint:
    fp_str = ''.join(map(byte_to_bits, fingerprint))
    return BitFingerprint(np.fromiter(map(int, fp_str), dtype=bool)[:bits_cnt])
