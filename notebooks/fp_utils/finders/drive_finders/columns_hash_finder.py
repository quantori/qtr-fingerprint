from typing import Iterable, Generator, List, Type, Optional
from abc import ABC, abstractmethod
import pandas as pd

from fp_utils.consts import PathType
from fp_utils.finders.finder import Finder
from fp_utils.finders.drive_finders.drive_finder import DriveFinder


class ColumnsHashFinder(DriveFinder, ABC):
    @property
    @abstractmethod
    def inner_finder_class(self) -> Type[Finder]:
        raise NotImplementedError

    @property
    @abstractmethod
    def columns(self) -> List[int]:
        raise NotImplementedError

    def __init__(self, df: pd.DataFrame, directory: PathType, unique_id: Optional[str] = None, *args, **kwargs) -> None:
        self.make_data_directory()
        hashes = pd.DataFrame(df.apply(self.__get_hash, axis=1), columns=['hash'])
        buckets = hashes.groupby('hash')
        self.inner_finders = dict()
        hash_values = buckets.apply(
            lambda bucket: self.__save_func(df.loc[bucket.index], bucket.iloc[0]['hash'], *args, **kwargs)).values
        self.hash_values = set(hash_values)

    def __save_func(self, bucket: pd.DataFrame, bucket_hash: int, *args, **kwargs) -> int:
        inner_finder = self.inner_finder_class(bucket, directory=self.data_directory, unique_id=str(bucket_hash),
                                               *args, **kwargs)
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
