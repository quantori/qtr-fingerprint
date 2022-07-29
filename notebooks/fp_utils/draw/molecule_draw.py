from typing import Tuple, Iterable
from rdkit import Chem
import rdkit.Chem.Draw
import matplotlib.pyplot as plt
from PIL import Image


class MoleculeDrawer:
    @staticmethod
    def draw_one(smiles: str) -> Image:
        return Chem.Draw.MolToImage(Chem.MolFromSmiles(smiles))

    @staticmethod
    def draw_many(smiles_list: Iterable[str], shape: Tuple[int, int] = (3, 3),
                  figsize: Tuple[float, float] = (16, 16)) -> None:
        images = [Chem.Draw.MolToImage(Chem.MolFromSmiles(smile)) for smile in smiles_list[:shape[0] * shape[1]]]
        _, axs = plt.subplots(shape[0], shape[1], figsize=figsize)
        axs = axs.flatten()
        for img, ax in zip(images, axs):
            ax.imshow(img)
