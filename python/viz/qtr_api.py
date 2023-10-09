from __future__ import annotations

import json
import os
from typing import TypedDict

import requests

API_HOST_ENV = 'QTR_API_HOST'
DEFAULT_API_HOST = "http://127.0.0.1"
API_PORT_ENV = 'QTR_API_PORT'
DEFAULT_API_PORT = 8080

CompoundSmiles = str  # type alias


class MoleculeRecord(TypedDict):
    libraryId: str
    id: CompoundSmiles


class SearchResults(TypedDict):
    isFinished: bool
    molecules: list[MoleculeRecord]


class QtrFingerprintApi:
    def __init__(self, host: str = DEFAULT_API_HOST, port: int = DEFAULT_API_PORT):
        self.host = host
        self.port = port

    def query_similar_compounds(self, compound: CompoundSmiles, limit: int) -> list[CompoundSmiles]:
        response = self._make_request(compound, limit)
        results = self._parse_response(response)
        return results

    def _make_request(self, target: CompoundSmiles, limit: int) -> SearchResults:
        url = f"{self.host}:{self.port}/fastquery"

        headers = {
            'Content-Type': 'application/json'
        }

        params = {
            'limit': limit,
        }

        payload = json.dumps(
            {
                "smiles": target,
                # not used yet
                "PUBCHEM_COMPONENT_COUNT": {
                    "min": 3,
                    "max": 3.1
                }
            }
        )
        response = requests.get(url=url, headers=headers, params=params, data=payload)
        response.raise_for_status()
        results = response.json()
        if results == "-1":
            raise APIError("Invalid response from server")
        return results

    def _parse_response(self, response: SearchResults) -> list[CompoundSmiles]:
        molecules = response['molecules']
        compounds = [mol['id'] for mol in molecules]
        return compounds

    @classmethod
    def from_env(cls) -> QtrFingerprintApi:
        host = os.getenv(API_HOST_ENV, DEFAULT_API_HOST)
        port = os.getenv(API_PORT_ENV, DEFAULT_API_PORT)
        return cls(host=host, port=port)


class APIError(Exception):
    pass


def test_query():
    api = QtrFingerprintApi.from_env()
    test_compound = 'c1cc(C(=O)O)c(OC(=O)C)cc1'
    results = api.query_similar_compounds(test_compound, 10)
    for result in results:
        print(result)


if __name__ == "__main__":
    test_query()
