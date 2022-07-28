import pandas as pd
from typing import Iterable
from fp_utils.finders.ram_finders.ram_finder import RamFinder
from fp_utils.settings import is_sub_fingerprint


class BFRamFinder(RamFinder):
    def __init__(self, df: pd.DataFrame, *args, **kwargs):
        self.df = df

    def find_all(self, fingerprint: pd.Series) -> Iterable[str]:
        for i, fp in enumerate(self.df.values):
            if is_sub_fingerprint(fingerprint, fp):
                yield self.df.index[i]
