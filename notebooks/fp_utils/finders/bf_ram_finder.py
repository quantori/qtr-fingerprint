import pandas as pd
import numpy as np
from typing import Iterable
from fp_utils.finders.ram_finder import RamFinder


class BFRamFinder(RamFinder):
    def __init__(self, df: pd.DataFrame):
        self.df = df

    def _find(self, fingerprint: pd.Series) -> Iterable[str]:
        for i, fp in enumerate(self.df.values):
            if np.all(fingerprint.values <= fp):
                yield self.df.index[i]
