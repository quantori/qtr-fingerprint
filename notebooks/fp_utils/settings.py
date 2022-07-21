from fp_utils import consts

from pandarallel import pandarallel
from fp_utils.consts import FingerprintType, FingerprintsArrayType
import pandas as pd
from typing import Iterable


def set_progress_bar(value: bool) -> None:
    pandarallel.initialize(progress_bar=value)


def fingerprint_to_series(fingerprint: FingerprintType) -> pd.Series:
    if type(fingerprint) == pd.DataFrame:
        fingerprint = fingerprint.iloc[0]
    return fingerprint


def fingerprint_array_to_iterable_of_series(fingerprints: FingerprintsArrayType) -> Iterable[FingerprintType]:
    if type(fingerprints) == pd.DataFrame:
        fingerprints = [fingerprints.iloc[i] for i in range(len(fingerprints))]
    return fingerprints


def init_fp_utils() -> None:
    set_progress_bar(consts.pandarallel_progress_bar)
