import pandas as pd
import numpy as np
from typing import Iterable
from pathlib import Path
from fp_utils.finders.drive_finder import DriveFinder
from fp_utils.consts import PathType


class BFDriveFinder(DriveFinder):
    def __init__(self, df: pd.DataFrame, directory: PathType, finder_name: str = "BFDriveFinder", *args, **kwargs):
        self.finder_path = Path(directory) / (finder_name + self.file_extension)
        self._pack(df, self.finder_path)

    def find_all(self, fingerprint: pd.Series) -> Iterable[str]:
        df = self._unpack(self.finder_path)
        for i, fp in enumerate(df.values):
            if np.all(fingerprint.values <= fp):
                yield df.index[i]
