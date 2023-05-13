import requests
import json
from enum import Enum
from abc import ABC, abstractmethod
from typing import List, Dict
from ipaddress import ip_address

from pathlib import Path
import sys

packages_dir = Path(__file__).absolute().parent.parent / 'packages'
assert packages_dir.is_dir()
sys.path.append(str(packages_dir))

from fp_utils import CatchTime


class SubstructureServer(ABC):
    @abstractmethod
    def ask(self, smiles: str, limit: int):
        raise NotImplementedError()


class PubchemServer(SubstructureServer):
    def ask(self, smiles: str, limit: int) -> List[int]:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsubstructure/smiles/{smiles}/cids/json"
        params = {"MaxRecords": limit}
        response = requests.post(url, params=params)
        assert response.ok, f"{response.text}"
        return response.json()["IdentifierList"]["CID"]


class QtrServer(SubstructureServer):

    def __init__(self, ip: ip_address, port: int):
        self.ip = ip
        self.port = port

    def ask(self, smiles: str, limit: int):
        url = f"http://{self.ip}:{self.port}/fastquery"
        data = {"smiles": smiles}
        params = {"limit": limit}
        headers = {"Content-type": "application/json"}
        response = requests.get(url, json=data, params=params, headers=headers)

        assert response.ok, f"Something wrong with server: {response.text}"
        json = response.json()
        return json["molecules"]


local_qtr_server = QtrServer(ip_address("0.0.0.0"), 8080)
remote_qtr_server = QtrServer(ip_address("18.189.162.73"), 8080)
pubchem_server = PubchemServer()


@CatchTime("query")
def do_query(server, smiles, limit):
    lst = server.ask(smiles, limit)
    print(f"For query: {smiles} found {len(lst)} answers")


@CatchTime("Total time")
def estimate_time(server, queries, limit):
    errors = 0
    successes = 0
    for smiles in queries:
        try:
            print(f"Start search for: {smiles}")
            do_query(server, smiles, limit)
            successes += 1
        except:
            print("Error occurred")
            errors += 1
    print(f"errors: {errors}, successes: {successes}")


def load_queries(filename):
    queries_path = Path(__file__).absolute().parent.parent.parent / "data" / filename
    queries = []
    with queries_path.open("r") as f:
        for line in f.readlines():
            smiles, _ = line.split()
            queries.append(smiles)
    return queries


if __name__ == "__main__":
    queries = load_queries("queries3488_sample_100.txt")
    estimate_time(pubchem_server, queries, 10000)

# pubchem
# errors: 12, successes: 88
# 278.009s -- Total time
# errors: 12, successes: 88
# 273.448s -- Total time


# qtr search
# errors: 0, successes: 100
# 213.895s -- Total time
# errors: 0, successes: 100
# 184.636s -- Total time
# errors: 0, successes: 100
# 195.417s -- Total time