from __future__ import annotations

from fp_utils.catch_time import CatchTime

from typing import List, Dict, Optional, Tuple
import numpy as np
import matplotlib.pyplot as plt


class SpeedTestStat:
    @property
    def measured_parameters(self) -> List[str]:
        return ['min', 'max', 'mean', 'median']

    def __init__(self, measurements: Optional[Dict[type, List[CatchTime]]] = None):
        self.measurements: Dict[type, List[CatchTime]] = measurements or dict()

    def __get_times(self, finder_class: type) -> List[float]:
        return list(map(lambda x: x.time, self.measurements[finder_class]))

    def __collect_stat(self) -> Dict[type, Dict[str, CatchTime]]:
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
        assert not self.measurements or set(self.measurements.keys()) == set(other.measurements.keys())
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
