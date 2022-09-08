import sys
from typing import Generator

from substructure_finder.db_filesystem import DbFilesystem
from substructure_finder.fingerprint import BitFingerprint
from substructure_finder import consts
from substructure_finder.splitter_tree import SplitterTree
from substructure_finder.search_engines.bucket_search_engine import BucketSearchEngine


class SearchEngine:
    def __init__(self, db_filesystem: DbFilesystem, db_name: str) -> None:
        self.splitter_tree = SplitterTree.load(db_filesystem.tree_path(db_name))
        self.pickle_paths = db_filesystem.pickle_paths(db_name)

    def search(self, fingerprint: BitFingerprint) -> Generator[str, None, None]:
        assert len(fingerprint) == consts.fingerprint_size_in_bits
        buckets = list(self.splitter_tree.get_buckets(fingerprint))
        print(f'Buckets number: {len(buckets)}', file=sys.stderr)
        for bucket in buckets:
            yield from BucketSearchEngine.search_in_file(fingerprint, self.pickle_paths[bucket])
