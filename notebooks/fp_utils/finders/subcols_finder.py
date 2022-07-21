from __future__ import annotations

import pandas as pd
import numpy as np
from typing import Iterable, Type, List
from abc import ABC, abstractmethod

from fp_utils.finders.finder import Finder


class SubColsFinder(Finder, ABC):
    COLUMNS = None

    @property
    def columns(self):
        if self.COLUMNS is None:
            raise NotImplementedError
        return self.COLUMNS

    def __init__(self, df: pd.DataFrame, base_finder: Type[Finder], *args, **kwargs) -> None:
        self.df = df
        self.sub_df = df.loc[df.index, df.columns.isin(self.columns)]
        self.sub_df = self.sub_df[self.columns]
        self.base_finder = base_finder(self.sub_df, *args, **kwargs)

    def _find(self, fingerprint: pd.Series) -> Iterable[str]:
        for answer in self.base_finder._find(fingerprint[self.columns]):
            if np.all(fingerprint.values <= self.df.loc[answer].values):
                yield answer
