from logging import getLogger
from typing import Optional

from indigo import Indigo, IndigoException, IndigoObject
from indigo.inchi import IndigoInchi
from indigo.renderer import IndigoRenderer
from rdkit import Chem

logger = getLogger(__name__)


def is_inchi(query) -> bool:
    return query.lower().startswith("inchi")


def render_result(indigo_renderer: IndigoRenderer, result: IndigoObject) -> bytes:
    try:
        image = indigo_renderer.renderToBuffer(result).tobytes()
    except IndigoException as e:
        logger.exception("Exception in rendering image")
        logger.exception(e)
        return b""
    return image


def highlight_key_in_query(indigo_instance: Indigo, result: IndigoObject, key: Optional[str]) -> IndigoObject:
    try:
        query_molecule = indigo_instance.loadQueryMolecule(key)
        query_molecule.aromatize()

        matcher = indigo_instance.substructureMatcher(result)
        match = matcher.match(query_molecule)

        if match:
            result = match.highlightedTarget()

    except IndigoException:
        logger.exception("Exception in loading query molecule, %s", key)

    return result


def render_image(query: str, key: Optional[str] = None) -> bytes:
    """Render image with indigo library.

    query: str - InChi or Smiles string to render image
    key: Optional[str]  - InChi or Smiles string to highlight in the result image
    """
    query_is_inchi = is_inchi(query)
    # Test the input query with RDKit
    mol = Chem.MolFromInchi(query) if query_is_inchi else Chem.MolFromSmiles(query)
    if not mol:
        return b""

    indigo_instance = Indigo()
    indigo_inchi_instance = IndigoInchi(indigo_instance)
    indigo_renderer = IndigoRenderer(indigo_instance)

    indigo_instance.setOption("render-output-format", "svg")
    indigo_instance.setOption("render-coloring", True)
    indigo_instance.setOption("render-implicit-hydrogens-visible", False)

    try:
        if query_is_inchi:
            molecule = indigo_inchi_instance.loadMolecule(query)
        else:
            molecule = indigo_instance.loadMolecule(query)
    except IndigoException as e:
        logger.exception("Exception in loading molecule, %s", query)
        logger.exception(e)
        return b""

    molecule.dearomatize()

    mol = None
    if key:
        # Test the input key with RDKit
        mol = Chem.MolFromInchi(key) if is_inchi(key) else Chem.MolFromSmiles(key)

    if mol:
        if is_inchi(key):
            key = indigo_inchi_instance.loadMolecule(key).smiles()

        molecule = highlight_key_in_query(indigo_instance, molecule, key)

    return render_result(indigo_renderer, molecule)
