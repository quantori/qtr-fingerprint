from __future__ import annotations

from typing import NewType
import numpy as np

BitFingerprint = NewType('Fingerprint', np.ndarray)
ByteFingerprint = NewType('ByteFingerprint', np.ndarray)
