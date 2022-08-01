from typing import Mapping, List, Iterable
from sklearn.neighbors import BallTree

from substrucure_finder.consts import Fingerprint
from substrucure_finder import utils
from substrucure_finder import consts


class BucketSearchEngine:
    def __init__(self, fingerprints: Mapping[str, Fingerprint], columns: List[int]) -> None:
        assert all(map(lambda x: len(x) == consts.fingerprint_size, fingerprints.values()))
        self.columns = columns
        self.id_to_smiles = dict()
        self.bucket_fingerprints = list()
        for i, (smiles, fingerprint) in enumerate(fingerprints.items()):
            self.id_to_smiles[i] = smiles
            self.bucket_fingerprints.append(fingerprint)
        sub_fingerprints = [utils.take_columns_from_fingerprint(self.columns, fp) for fp in self.bucket_fingerprints]
        self.ball_tree = BallTree(sub_fingerprints, leaf_size=2, metric='russellrao')

    def search(self, fingerprint: Fingerprint) -> Iterable[str]:
        sub_fingerprint = utils.take_columns_from_fingerprint(self.columns, fingerprint)
        radius = utils.russelrao_radius(sub_fingerprint)
        ball_tree_answers = self.ball_tree.query_radius([sub_fingerprint], radius)[0]
        filtered_answers = [i for i in ball_tree_answers if
                            utils.is_sub_fingerprint(fingerprint, self.bucket_fingerprints[i])]
        print(filtered_answers)
        print(ball_tree_answers)
        return [self.id_to_smiles[i] for i in filtered_answers]
