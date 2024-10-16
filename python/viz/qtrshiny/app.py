import hashlib
import logging
import time
from itertools import chain
from pathlib import Path
from typing import Callable, TypeVar

from rdkit import Chem
from rdkit.Chem import Draw
from shiny import App, Inputs, Outputs, Session, reactive, render, ui

from qtr_api import CompoundSmiles, QtrFingerprintApi

logger = logging.getLogger(__name__)

ROW_SIZE = 5
IMG_STYLE = "width: 250px;"
MAX_RESULTS = 200

DATABASE_SIZE = "240k"

app_ui = ui.page_fluid(

    ui.row(
        ui.column(
            3, ui.div(
                ui.output_text("database_size"),
                ui.input_slider("n", "Max number of compounds", min=1, max=MAX_RESULTS, value=1),
                ui.input_text("compound", "Compound:", 'c1cc(C(=O)O)c(OC(=O)C)cc1'),
                ui.input_action_button("search", "Substructure search", class_="btn-success"),
                ui.p("Target compound:"),
                ui.output_ui("target_image"),
            )
        ),
        ui.column(
            9,
            ui.output_text("result_count"),
            ui.output_text("search_time"),
            ui.output_text("visualisation_time"),
            ui.output_ui("compounds_grid"),
        ),
    ),
)


# noinspection PyUnusedLocal
def server(inputs: Inputs, outputs: Outputs, session: Session):
    @outputs
    @render.ui
    @reactive.event(inputs.compound)
    def target_image():
        return create_image_from_smiles(inputs.compound())

    @outputs
    @render.ui
    def compounds_grid():
        images, _ = result_image_files()
        # spits in chunks
        row_size = ROW_SIZE
        images = [images[i: i + row_size] for i in range(0, len(images), row_size)]
        images = [ui.div(row) for row in images]
        images = ui.div(images)
        return images

    @outputs
    @render.text
    def search_time():
        _, execution_time = compound_search()
        return f"Substructure search time: {execution_time:.4f} sec"

    @outputs
    @render.text
    def visualisation_time():
        _, render_time = result_image_files()
        return f"Molecules render time: {render_time:.4f} sec"

    @outputs
    @render.text
    def result_count():
        results, _ = compound_search()
        return f"Found {len(results)} compounds"

    @reactive.Calc
    @reactive.event(inputs.search, ignore_none=False, ignore_init=True)
    def compound_search() -> tuple[list[CompoundSmiles], float]:
        n_results: int = inputs.n()
        target_compound: CompoundSmiles = inputs.compound()

        logger.info(f"Performing substructure search for {target_compound},  max {n_results}")
        results, elapsed = timed(api.query_similar_compounds, target_compound, limit=n_results)
        logger.info(f"Search took {elapsed} seconds, returned {len(results)} results")

        return results, elapsed

    @reactive.Calc
    def result_image_files():
        results, _ = compound_search()
        # takes the input compound during last search
        # doesn't redraw when changing compounds
        with reactive.isolate():
            highlight = inputs.compound()
        logger.info(f"Creating {len(results)} images from search results")
        images, render_time = timed(
            lambda res: [create_image_from_smiles(compound, highlight) for compound in res],
            results
        )
        logger.info(f"Rendering result images took {render_time} seconds")
        return images, render_time

    def create_image_from_smiles(compound: CompoundSmiles, highlight: CompoundSmiles = None):
        filename = generate_filename(compound)
        filepath = www_dir / filename
        logger.debug(f"Turning compound {compound} into image at {filepath}")
        generate_image(filepath, compound, highlight)
        shown_image = ui.img(name=compound, src=filename, style=IMG_STYLE, border="1")
        return shown_image

    @outputs
    @render.text
    def database_size() -> str:
        return f"Database size is {DATABASE_SIZE}"


def generate_filename(compound: CompoundSmiles) -> str:
    # SMILES have symbols that lead to invalid filenames
    # we need them to be unique
    # BLAKE2 has enough collision resistance and fast enough
    compound_hash = hashlib.blake2b(compound.encode('utf-8')).hexdigest()
    return f'{compound_hash}.png'


def generate_image(file: Path, compound: CompoundSmiles, highlight: CompoundSmiles = None):
    logger.debug(f"Saving {compound} to {file}")
    mol = Chem.MolFromSmiles(compound)
    if highlight:
        part = Chem.MolFromSmiles(highlight)
        matches = mol.GetSubstructMatches(part)
        if matches:
            logger.debug(f"Found {len(matches)} matches in {compound}")
            match = chain(*matches)
        else:
            logger.warning(f"Unable to find {highlight} in {compound}")
            match = []
    else:
        match = []

    Draw.MolToFile(mol, file, highlightAtoms=match)


T = TypeVar('T')


def timed(func: Callable[..., T], *args, **kwargs) -> tuple[T, float]:
    """
    Returns result of function with its execution time
    """
    start_time = time.perf_counter()
    results = func(*args, **kwargs)
    finish_time = time.perf_counter()
    elapsed = finish_time - start_time
    return results, elapsed


logging.basicConfig(level=logging.INFO)

www_dir = Path(__file__).parent.parent / "www"
www_dir.mkdir(parents=True, exist_ok=True)
api = QtrFingerprintApi.from_env()
app = App(app_ui, server, static_assets=www_dir)
