from abc import ABC

import pandas as pd
from typing import Iterable, Optional
from fp_utils.finders.drive_finders.drive_finder import DriveFinder
from fp_utils.consts import PathType
from fp_utils.settings import is_sub_fingerprint


class BFDriveFinder(DriveFinder, ABC):
    def __init__(self, df: pd.DataFrame, directory: PathType, unique_id: Optional[str] = None, *args, **kwargs):
        self.make_data_directory()
        self._pack(df, self.data_directory)

    def find_all(self, fingerprint: pd.Series) -> Iterable[str]:
        df = self._unpack(self.data_directory)
        for i, fp in enumerate(df.values):
            if is_sub_fingerprint(fingerprint, fp):
                yield df.index[i]
