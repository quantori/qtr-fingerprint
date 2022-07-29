from __future__ import annotations

import pandas as pd
from typing import Iterable, Type, Optional
from abc import ABC, abstractmethod

from fp_utils.finders.drive_finders.drive_finder import DriveFinder
from fp_utils.finders.finder import Finder
from fp_utils.consts import PathType
from fp_utils.settings import is_sub_fingerprint


class SubColsFinder(DriveFinder, ABC):

    @property
    @abstractmethod
    def inner_finder_class(self) -> Type[Finder]:
        raise NotImplementedError

    @property
    @abstractmethod
    def columns(self) -> Iterable[int]:
        raise NotImplementedError

    def __init__(self, df: pd.DataFrame, directory: PathType, unique_id: Optional[str] = None, *args, **kwargs) -> None:
        self.make_data_directory()
        self.packer.pack(df, self.index_file)
        sub_df = df.loc[df.index, df.columns.isin(self.columns)]
        sub_df = sub_df[self.columns]
        self.inner_finder = self.inner_finder_class(sub_df, self.data_directory, unique_id=None, *args, **kwargs)

    def find_all(self, fingerprint: pd.Series) -> Iterable[str]:
        df = self.packer.unpack(self.index_file)
        for answer in self.inner_finder.find_all(fingerprint[self.columns]):
            if is_sub_fingerprint(fingerprint.values, df.loc[answer]):
                yield answer
