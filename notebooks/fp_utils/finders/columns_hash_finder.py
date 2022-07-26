from typing import Iterable, Generator, List, Type
from abc import ABC, abstractmethod
import pandas as pd
from pathlib import Path

from fp_utils.consts import PathType
from fp_utils.finders.finder import Finder


class ColumnsHashFinder(Finder, ABC):
    @property
    @abstractmethod
    def inner_finder_class(self) -> Type[Finder]:
        raise NotImplementedError

    @property
    @abstractmethod
    def columns(self) -> List[int]:
        raise NotImplementedError

    def __init__(self, df: pd.DataFrame, directory: PathType, finder_name: str, *args, **kwargs) -> None:
        self.finder_path = Path(directory) / finder_name
        self.finder_path.mkdir(parents=True, exist_ok=True)
        hashes = pd.DataFrame(df.apply(self.__get_hash, axis=1), columns=['hash'])
        buckets = hashes.groupby('hash')
        self.inner_finders = dict()
        hash_values = buckets.apply(
            lambda bucket: self.__save_func(df.loc[bucket.index], bucket.iloc[0]['hash'], *args, **kwargs)).values
        self.hash_values = set(hash_values)

    def __save_func(self, bucket: pd.DataFrame, bucket_hash: int, *args, **kwargs) -> int:
        inner_finder = self.inner_finder_class(bucket, directory=self.finder_path, finder_name=str(bucket_hash), *args,
                                               **kwargs)
        self.inner_finders[bucket_hash] = inner_finder
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

    def find_all(self, fingerprint: pd.Series) -> Iterable[str]:
        query_mask = self.__get_hash(fingerprint)
        for bucket_mask in self.__get_meta_masks(query_mask):
            if bucket_mask not in self.hash_values:
                continue
            inner_finder = self.inner_finders[bucket_mask]
            yield from inner_finder.find_all(fingerprint)
