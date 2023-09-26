from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Draw
from shiny import App, Inputs, Outputs, Session, reactive, render, ui

ROW_SIZE = 5
IMG_STYLE = "width: 250px;"

CompoundSmiles = str  # type alias

app_ui = ui.page_fluid(

    ui.row(
        ui.column(
            3, ui.div(
                ui.p("Database size is 112kk"),
                ui.input_slider("n", "Max number of compounds", min=0, max=27, value=1),
                ui.input_text("compound", "Compound:", 'c1cc(C(=O)O)c(OC(=O)C)cc1'),
                ui.input_action_button("search", "Substructure search", class_="btn-success"),
                ui.p("Target compound:"),
                ui.output_ui("target_image"),
            )
        ),
        ui.column(
            9,
            ui.output_ui("compounds_grid"),
        ),
    ),
)


def server(inputs: Inputs, outputs: Outputs, session: Session):
    @outputs
    @render.ui
    @reactive.event(inputs.compound)
    def target_image():
        return create_image_from_smiles(inputs.compound())

    @outputs
    @render.ui
    @reactive.event(inputs.search, ignore_none=False)
    def compounds_grid():
        n_results: int = inputs.n()
        target_compound: CompoundSmiles = inputs.compound()

        similar_compounds = find_similar_compounds(target_compound, limit=n_results)
        result_grid = create_image_grid(similar_compounds, ROW_SIZE)
        return result_grid

    def create_image_grid(results: list[str], row_size: int):
        images = [create_image_from_smiles(compound) for compound in results]
        # spits in chunks
        images = [images[i: i + row_size] for i in range(0, len(results), row_size)]
        images = [ui.div(row) for row in images]
        images = ui.div(images)
        return images

    def create_image_from_smiles(compound: CompoundSmiles):
        filename = f"{compound}.png"
        filepath = www_dir / filename
        generate_image(compound, file=filepath)
        shown_image = ui.img(name=compound, src=filename, style=IMG_STYLE, border="1")
        return shown_image


def find_similar_compounds(target: CompoundSmiles, limit: int = None) -> list[CompoundSmiles]:
    # stub
    if limit is None:
        limit = 13
    return [target] * limit


def generate_image(smiles: CompoundSmiles, file: Path):
    mol = Chem.MolFromSmarts(smiles)
    Draw.MolToFile(mol, file)


www_dir = Path(__file__).parent / "www"
www_dir.mkdir(parents=True, exist_ok=True)
app = App(app_ui, server, static_assets=www_dir)
