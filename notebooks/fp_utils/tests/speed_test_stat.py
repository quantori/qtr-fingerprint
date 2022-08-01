from __future__ import annotations

from fp_utils.catch_time import CatchTime
from fp_utils.finders.finder import Finder
from typing import List, Dict, Optional, Tuple, Iterable
import numpy as np
import matplotlib.pyplot as plt


class SpeedTestStat:
    @property
    def measured_parameters(self) -> List[str]:
        return ['min', 'max', 'mean', 'median']

    def __init__(self, measurements: Optional[Dict[Finder, List[float]]] = None):
        self.measurements: Dict[Finder, List[CatchTime]] = measurements or dict()

    def __collect_stat(self) -> Dict[Finder, Dict[str, CatchTime]]:
        stat = dict()
        for finder, measurements in self.measurements.items():
            stat[finder] = dict()
            for parameter in self.measured_parameters:
                value = getattr(np, parameter)(measurements)
                stat[finder][parameter] = value
        return stat

    def get_name(self, finder: Finder) -> str:
        name = finder.__class__.__name__
        cnt = 0
        for other_finder in self.measurements.keys():
            cnt += int(name == other_finder.__class__.__name__)
        assert cnt >= 1
        return name if cnt == 1 else f'{name}_{id(finder)}'

    def __str__(self):
        stat = self.__collect_stat()
        res = []
        for finder in stat:
            res += [f'{self.get_name(finder)}:']
            for parameter in self.measured_parameters:
                res += [f'\t{stat[finder][parameter]:.3f} -- {parameter}']
        return '\n'.join(res)

    def __repr__(self):
        return str(self)

    def __iadd__(self, other: SpeedTestStat) -> SpeedTestStat:
        for finder, measurements in other.measurements.items():
            self.measurements.setdefault(finder, [])
            self.measurements[finder] += measurements
        return self

    def __add__(self, other: SpeedTestStat) -> SpeedTestStat:
        res = SpeedTestStat()
        res += self
        res += other
        return res

    def as_plot(self, figsize: Tuple[float, float] = (16, 8)) -> None:
        assert self.measurements
        _, ax = plt.subplots(figsize=figsize)
        for finder, measurements in self.measurements.items():
            ax.plot(measurements, label=self.get_name(finder))
        ax.legend()
        plt.show()

    def as_boxplot(self, figsize: Tuple[float, float] = (8, 8)) -> None:
        assert self.measurements
        _, ax = plt.subplots(figsize=figsize)
        data = [self.measurements[key] for key in self.measurements]
        labels = [self.get_name(finder) for finder in self.measurements]
        ax.boxplot(data, vert=True, patch_artist=True, labels=labels)
        plt.show()

    def drop(self, finders: Iterable[Finder]) -> SpeedTestStat:
        good_finders = set(self.measurements.keys()) - set(finders)
        return self.take(good_finders)

    def take(self, finders: Iterable[Finder]) -> SpeedTestStat:
        new_stat = SpeedTestStat()
        for finder in finders:
            new_stat.measurements[finder] = self.measurements[finder]
        return new_stat
