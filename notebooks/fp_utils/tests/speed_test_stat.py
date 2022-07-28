from __future__ import annotations

from fp_utils.catch_time import CatchTime
from fp_utils.finders.finder import Finder
from typing import List, Dict, Optional, Tuple, Type, Union
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy


class SpeedTestStat:
    @property
    def measured_parameters(self) -> List[str]:
        return ['min', 'max', 'mean', 'median']

    def __init__(self, measurements: Optional[Dict[str, List[CatchTime]]] = None):
        self.measurements: Dict[str, List[CatchTime]] = measurements or dict()

    def __get_times(self, finder_name: str) -> List[float]:
        return list(map(lambda x: x.time, self.measurements[finder_name]))

    def __collect_stat(self) -> Dict[str, Dict[str, CatchTime]]:
        stat = dict()
        for finder_name in self.measurements:
            times = np.array(self.__get_times(finder_name))
            stat[finder_name] = dict()
            for parameter in self.measured_parameters:
                value = getattr(np, parameter)(times)
                stat[finder_name][parameter] = value
        return stat

    def __str__(self):
        stat = self.__collect_stat()
        res = []
        for finder_name in stat:
            res += [f'{finder_name}:']
            for parameter in self.measured_parameters:
                res += [f'\t{stat[finder_name][parameter]:.3f} -- {parameter}']
        return '\n'.join(res)

    def __repr__(self):
        return str(self)

    def __iadd__(self, other: SpeedTestStat) -> SpeedTestStat:
        for key in other.measurements:
            self.measurements.setdefault(key, [])
            self.measurements[key] += other.measurements[key]
        return self

    def __add__(self, other: SpeedTestStat) -> SpeedTestStat:
        res = SpeedTestStat()
        res += self
        res += other
        return res

    def as_plot(self, figsize: Tuple[float, float] = (16, 8)) -> None:
        assert self.measurements
        _, ax = plt.subplots(figsize=figsize)
        for key in self.measurements:
            ax.plot(self.__get_times(key), label=key.__name__)
        ax.legend()
        plt.show()

    def as_boxplot(self, figsize: Tuple[float, float] = (8, 8)) -> None:
        assert self.measurements
        _, ax = plt.subplots(figsize=figsize)
        data = [self.__get_times(key) for key in self.measurements]
        labels = [key.__name__ for key in self.measurements]
        ax.boxplot(data, vert=True, patch_artist=True, labels=labels)
        plt.show()

    def drop(self, finders: List[Finder]) -> SpeedTestStat:
        new_stat = deepcopy(self)
        for finder in finders:
            new_stat.measurements.pop(finder.__name__, None)
        return new_stat

    def take(self, finders: List[Finder]) -> SpeedTestStat:
        new_stat = SpeedTestStat()
        for finder in finders:
            new_stat.measurements[finder.__name__] = self.measurements[finder.__name__]
        return new_stat
