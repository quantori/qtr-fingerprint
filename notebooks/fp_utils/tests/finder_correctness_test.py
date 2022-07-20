import numpy as np
import pandas as pd
import random

from typing import Iterable

from fp_utils.logger import Logger
from fp_utils.finders.finder import Finder
from fp_utils.consts import FingerprintType, FingerprintsArrayType
from fp_utils import settings


class FinderCorrectnessTester:
    def __init__(self, correct_finder: Finder, finders: Iterable[Finder]) -> None:
        self.correct_finder = correct_finder
        self.finders = list(finders)

    def test_one(self, fingerprint: FingerprintType, verbose: bool = False) -> np.ndarray:
        fingerprint = settings.fingerprint_to_series(fingerprint)
        correct_results = set(self.correct_finder.find(fingerprint))
        answers = np.array([True] * len(self.finders))
        for i, finder in enumerate(self.finders):
            results = set(finder.find(fingerprint))
            cmp = results == correct_results
            answers[i] &= cmp
            Logger.log(f'{"OK" if cmp else "WA"} -- {finder.__class__.__name__}', verbose)
        return answers

    def test_all(self, fingerprints: FingerprintsArrayType, verbose: bool = False) -> np.ndarray:
        fingerprints = settings.fingerprint_array_to_iterable_of_series(fingerprints)
        answers = np.array([True] * len(self.finders))
        Logger.log('-----', verbose)
        for i, fingerprint in enumerate(fingerprints):
            Logger.log(f'Test #{i:03d}', verbose)
            answers &= self.test_one(fingerprint, verbose)
            Logger.log('-----', verbose)
        return answers

    def test_random(self, fingerprints: FingerprintsArrayType, test_count: int, verbose: bool = False) -> np.ndarray:
        if isinstance(fingerprints, pd.DataFrame):
            return self.test_all(fingerprints.sample(test_count), verbose)
        else:
            return self.test_all(random.sample(list(fingerprints), test_count), verbose)
