import pandas as pd
from typing import Iterable
from fp_utils.finders.ram_finders.ram_finder import RamFinder
from fp_utils.settings import is_sub_fingerprint


class BFPRamFinder(RamFinder):
    def __init__(self, df: pd.DataFrame, *args, **kwargs):
        self.df = df

    def find_all(self, fingerprint: pd.Series) -> Iterable[str]:
        d = self.df.parallel_apply(lambda row: is_sub_fingerprint(fingerprint.values, row), axis=1)
        return d[d].index
