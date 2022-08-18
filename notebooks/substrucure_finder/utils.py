from typing import List
import numpy as np
from pathlib import Path

from substrucure_finder.fingerprint import Fingerprint
from substrucure_finder import consts


def fingerprint_to_str(fingerprint: Fingerprint) -> str:
    return ''.join(map(str, fingerprint))


def russelrao_radius(fingerprint: Fingerprint) -> float:
    return 1 - sum(fingerprint) / len(fingerprint) + consts.radius_eps


def take_columns_from_fingerprint(columns: List[int], fingerprint: Fingerprint) -> Fingerprint:
    return Fingerprint(np.array([fingerprint[i] for i in columns]))


def is_sub_fingerprint(sub_fingerprint: Fingerprint, meta_fingerprint: Fingerprint) -> bool:
    assert len(sub_fingerprint) == len(meta_fingerprint)
    return np.all(sub_fingerprint <= meta_fingerprint)


def bucket_path(data_path: Path, bucket: int) -> Path:
    return data_path / 'buckets' / (str(bucket) + '.pickle')


def raw_bucket_path(data_path: Path, bucket: int) -> Path:
    return data_path / 'raw_buckets' / str(bucket)


def splitter_tree_path(data_path: Path) -> Path:
    return data_path / 'tree'


def columns_path(data_path: Path, bucket: int) -> Path:
    return data_path / 'raw_buckets' / (str(bucket) + '.col')


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