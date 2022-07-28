import pickle
from pathlib import Path

from substrucure_finder.bucket_search_engine import BucketSearchEngine


class Dumper:

    def __init__(self, file: Path) -> None:
        self.file = file

    def __enter__(self):
        self.stream = self.file.open('wb')

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.stream.close()

    def bucket_search_engine(self, bucket_search_engine: BucketSearchEngine) -> None:
        pickle.dump(bucket_search_engine, self.stream, protocol=-1)
