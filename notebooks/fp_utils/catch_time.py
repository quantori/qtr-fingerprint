import time
from enum import IntEnum
from typing import Callable
import functools


class CatchTimeState(IntEnum):
    NOT_STARTED = 0
    IN_PROGRESS = 1
    FINISHED = 2


class CatchTime:
    def __init__(self, text='CatchTime'):
        self.text = text
        self.state = CatchTimeState.NOT_STARTED

    def __enter__(self):
        self.time = time.perf_counter()
        self.state = CatchTimeState.IN_PROGRESS
        return self

    def __exit__(self, *kwargs):
        self.time = time.perf_counter() - self.time
        self.state = CatchTimeState.FINISHED

    def __str__(self):
        if self.state == CatchTimeState.FINISHED:
            return f'{self.time:.3f} -- {self.text}'
        else:
            return f'{self.state.name} -- {self.text}'

    def __repr__(self):
        if self.state == CatchTimeState.FINISHED:
            return f'{self.__class__.__name__}({self.text}, {self.time:.3f} sec)'
        else:
            return f'{self.__class__.__name__}({self.text}, {self.state.name})'

    def __call__(self, function: Callable) -> Callable:
        @functools.wraps(function)
        def wrap_function(*args, **kwargs):
            with self:
                return_value = function(*args, **kwargs)
            print(self)
            return return_value

        return wrap_function
