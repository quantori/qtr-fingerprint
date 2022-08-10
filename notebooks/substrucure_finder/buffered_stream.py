from pathlib import Path
import sys


class BufferedStream:
    def __init__(self, file_path: Path):
        with file_path.open('rb') as stream:
            # print("start read buffer", file=sys.stderr)
            self.buffer = stream.read()
            # print("finish read buffer", file=sys.stderr)
        self.ptr = 0

    def read(self, n: int):
        l_ptr = self.ptr
        self.ptr = min(self.ptr + n, len(self.buffer))
        return self.buffer[l_ptr:self.ptr]

    def read_until(self, symbol: str):
        b = symbol.encode()
        l_ptr = self.ptr
        self.ptr = self.buffer.find(b, self.ptr)
        assert self.ptr != -1
        self.ptr += 1
        # l_ptr = self.ptr
        # while True:
        #     self.ptr += 1
        #     assert self.ptr < len(self.buffer)
        #     if self.buffer[self.ptr] == s:
        #         break
        # assert chr(self.buffer[self.ptr]) == symbol
        # self.ptr += 1
        return self.buffer[l_ptr:self.ptr - 1]
