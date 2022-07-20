import pandas as pd
import numpy as np
from typing import Iterable, Any, BinaryIO
from pathlib import Path
import pickle
from fp_utils.finders.drive_finder import DriveFinder
from fp_utils.consts import PathType


class BFDriveFinder(DriveFinder):
    @property
    def df_name(self):
        return "df.pickle"

    def __init__(self, df: pd.DataFrame, directory: PathType):
        self.df_path = Path(directory) / self.df_name
        self._pack(df, self.df_path)

    def _find(self, fingerprint: pd.Series) -> Iterable[str]:
        df = self._unpack(self.df_path)
        for i, fp in enumerate(df.values):
            if np.all(fingerprint.values <= fp):
                yield df.index[i]

    def _dump(self, obj: Any, file: BinaryIO) -> None:
        pickle.dump(obj, file, protocol=-1)

    def _load(self, file: BinaryIO) -> Any:
        return pickle.load(file)
