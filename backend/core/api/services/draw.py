import io

from rdkit import Chem
from rdkit.Chem.Draw import MolDrawOptions, rdMolDraw2D


def draw_molecule(mol, r, g, b, molSize=(450, 150), kekulize=True) -> io.StringIO:
    mc = Chem.MolFromSmiles(mol)
    if kekulize:
        try:
            Chem.Kekulize(mc)
        except Exception as e:
            mc = Chem.Mol(mol.ToBinary())

    if not mc.GetNumConformers():
        Chem.rdDepictor.Compute2DCoords(mc)

    optn = MolDrawOptions()
    optn.setBackgroundColour((r, g, b))
    drawer = rdMolDraw2D.MolDraw2DSVG(*molSize)
    drawer.SetDrawOptions(optn)
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    buffer = io.StringIO()
    buffer.write(drawer.GetDrawingText())
    buffer.seek(0)
    return buffer
