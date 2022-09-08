from typing import Union, Iterable
from pathlib import Path
import pandas as pd

pandarallel_progress_bar = True

FingerprintType = Union[pd.Series, pd.DataFrame]
FingerprintsArrayType = Union[Iterable[FingerprintType], pd.DataFrame]
PathType = Union[str, Path]
