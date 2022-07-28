from typing import Literal, NewType
import numpy as np

byteorder: Literal["little", "big"] = "little"
uint64_minus_one = int.from_bytes(bytes([0xff] * 8), byteorder, signed=False)

Fingerprint = NewType('Fingerprint', np.ndarray)

fingerprint_size = 2584  # TODO: set correct value
fingerprint_size_in_bytes = (fingerprint_size + 7) / 8
columns_count = 64
radius_eps = 1e-5
