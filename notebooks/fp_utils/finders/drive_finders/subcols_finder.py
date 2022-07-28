from __future__ import annotations

import pandas as pd
from typing import Iterable, Type
from abc import ABC, abstractmethod
from pathlib import Path

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

    @property
    def df_path(self) -> Path:
        return self.finder_path / ('df' + self.file_extension)

    def __init__(self, df: pd.DataFrame, directory: PathType, finder_name: str, *args, **kwargs) -> None:
        self.finder_path = Path(directory) / finder_name
        self.finder_path.mkdir(parents=True, exist_ok=True)
        self._pack(df, self.df_path)
        sub_df = df.loc[df.index, df.columns.isin(self.columns)]
        sub_df = sub_df[self.columns]
        self.inner_finder = self.inner_finder_class(sub_df, self.finder_path, self.inner_finder_class.__name__, *args,
                                                    **kwargs)

    def find_all(self, fingerprint: pd.Series) -> Iterable[str]:
        df = self._unpack(self.df_path)
        for answer in self.inner_finder.find_all(fingerprint[self.columns]):
            if is_sub_fingerprint(fingerprint, df.loc[answer]):
                yield answer
