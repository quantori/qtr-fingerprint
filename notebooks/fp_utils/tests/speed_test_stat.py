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

    def __init__(self, measurements: Optional[Dict[Type[Finder], List[CatchTime]]] = None):
        self.measurements: Dict[Type[Finder], List[CatchTime]] = measurements or dict()

    def __get_times(self, finder_class: Type[Finder]) -> List[float]:
        return list(map(lambda x: x.time, self.measurements[finder_class]))

    def __collect_stat(self) -> Dict[Type[Finder], Dict[str, CatchTime]]:
        stat = dict()
        for finder_class in self.measurements:
            times = np.array(self.__get_times(finder_class))
            stat[finder_class] = dict()
            for parameter in self.measured_parameters:
                value = getattr(np, parameter)(times)
                stat[finder_class][parameter] = value
        return stat

    def __str__(self):
        stat = self.__collect_stat()
        res = []
        for engine_class in stat:
            res += [f'{engine_class.__name__}:']
            for parameter in self.measured_parameters:
                res += [f'\t{stat[engine_class][parameter]:.3f} -- {parameter}']
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

    def drop(self, finders: List[Union[Finder, Type[Finder]]]) -> SpeedTestStat:
        new_stat = deepcopy(self)
        for finder in finders:
            if isinstance(finder, Finder):
                finder = finder.__class__
            assert issubclass(finder, Finder)
            new_stat.measurements.pop(finder, None)
        return new_stat

    def take(self, finders: List[Union[Finder, Type[Finder]]]) -> SpeedTestStat:
        new_stat = SpeedTestStat()
        for finder in finders:
            if isinstance(finder, Finder):
                finder = finder.__class__
            assert issubclass(finder, Finder)
            new_stat.measurements[finder] = self.measurements[finder]
        return new_stat
