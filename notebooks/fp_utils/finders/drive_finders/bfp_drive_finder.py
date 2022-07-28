import pandas as pd
from typing import Iterable
from pathlib import Path
from fp_utils.finders.drive_finders.drive_finder import DriveFinder
from fp_utils.consts import PathType
from fp_utils.settings import is_sub_fingerprint


class BFPDriveFinder(DriveFinder):
    def __init__(self, df: pd.DataFrame, directory: PathType, finder_name: str = "BFPDriveFinder", *args, **kwargs):
        self.finder_path = Path(directory) / (finder_name + self.file_extension)
        self._pack(df, self.finder_path)

    def find_all(self, fingerprint: pd.Series) -> Iterable[str]:
        df = self._unpack(self.finder_path)
        d = df.parallel_apply(lambda row: is_sub_fingerprint(fingerprint, row), axis=1)
        return d[d].index
