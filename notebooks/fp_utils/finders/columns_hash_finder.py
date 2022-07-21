import pickle
from typing import Any, BinaryIO, Iterable, Generator, List
from pathlib import Path
from abc import ABC, abstractmethod
import pandas as pd
from sklearn.neighbors import BallTree
from pandarallel import pandarallel

from fp_utils.consts import PathType
from fp_utils.finders.drive_finder import DriveFinder


class ColumnsHashFinder(DriveFinder, ABC):
    FILE_EXTENSION = '.pickle'
    R_EPS = 1e-5

    @property
    @abstractmethod
    def columns(self) -> List[int]:
        raise NotImplementedError

    def __init__(self, df: pd.DataFrame, directory: PathType) -> None:
        self.directory = Path(directory)
        hashes = pd.DataFrame(df.apply(self.__get_hash, axis=1), columns=['hash'])
        buckets = hashes.groupby('hash')
        hash_values = buckets.apply(lambda bucket: self.__save_func(bucket, df)).values
        self.hash_values = set(hash_values)

    def __get_bucket_path(self, bucket_hash: int) -> Path:
        return self.directory / (str(bucket_hash) + self.FILE_EXTENSION)

    def __save_func(self, bucket: pd.DataFrame, df: pd.DataFrame) -> int:
        index = bucket.index, BallTree(df.loc[bucket.index].values, leaf_size=2, metric='russellrao')
        bucket_hash = bucket.iloc[0]['hash']
        path = self.__get_bucket_path(bucket_hash)
        self._pack(index, path)
        return bucket_hash

    def __get_hash(self, fingerprint: pd.Series) -> int:
        h = 0
        for c in self.columns:
            h = (h << 1) + fingerprint[c]
        return h

    def __get_meta_masks(self, mask: int) -> Generator[int, None, None]:
        x = (1 << len(self.columns)) - 1
        t = mask ^ x
        while True:
            meta_mask = t ^ x
            yield meta_mask
            if t == 0:
                break
            t = (t - 1) & (mask ^ x)

    def _find(self, fingerprint: pd.Series) -> Iterable[str]:
        query_mask = self.__get_hash(fingerprint)
        r = 1 - sum(fingerprint) / len(fingerprint) + self.R_EPS

        for bucket_mask in self.__get_meta_masks(query_mask):
            if bucket_mask not in self.hash_values:
                continue

            path = self.__get_bucket_path(bucket_mask)
            bucket_index, tree = self._unpack(path)
            ind = tree.query_radius([fingerprint], r)
            for pos in ind[0]:
                yield bucket_index[pos]

    def _load(self, file: BinaryIO) -> Any:
        return pickle.load(file)

    def _dump(self, obj: Any, file: BinaryIO) -> None:
        pickle.dump(obj, file)
