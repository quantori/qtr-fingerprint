from typing import Literal

byteorder: Literal["little", "big"] = "little"
uint64_minus_one = int.from_bytes(bytes([0xff] * 8), byteorder, signed=False)

fingerprint_size_in_bits = 2584
fingerprint_size_in_bytes = (fingerprint_size_in_bits + 7) // 8
columns_count = 64
radius_eps = 1e-5
