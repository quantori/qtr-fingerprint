from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Draw
from shiny import App, Inputs, Outputs, Session, reactive, render, ui

ROW_SIZE = 5
IMG_STYLE = "width: 250px;"


app_ui = ui.page_fluid(

    ui.row(
        ui.column(
            3, ui.div(
                ui.p("Database size is 112kk"),
                ui.input_slider("n", "Max number of compounds", min=0, max=1000, value=1),
                ui.input_text("compound", "Compound:", 'c1cc(C(=O)O)c(OC(=O)C)cc1'),
                ui.input_action_button("search", "Substructure search", class_="btn-success"),
                ui.p("Target compount:"),
                ui.output_ui("target"),
            )
        ),
        ui.column(
            9,
            ui.output_ui("result"),
        ),
    ),
)


def square(x: ui.TagChild, n: int) -> ui.Tag:
    row = ui.div([x] * n)
    return ui.div([row] * n)


def server(input: Inputs, output: Outputs, session: Session):
    def save_smiles(smiles):
        mol = Chem.MolFromSmarts(smiles)
        file_name = smiles + ".png"
        Draw.MolToFile(mol, www_dir / file_name)
        return file_name

    def generate_smiles_img(smiles, **metadata):
        file_name = save_smiles(smiles)
        return ui.img(name=smiles, src=file_name, style=IMG_STYLE, **metadata)

    @output
    @render.ui
    @reactive.event(input.search, ignore_none=False)
    def target():
        return generate_smiles_img(input.compound())

    def found_smiles(smiles):
        return [smiles] * 13

    @output
    @render.ui
    @reactive.event(input.search, ignore_none=False)
    def result():
        smiless = found_smiles(input.compound())
        imgs = []
        row = []
        for smiles in smiless:
            row.append(generate_smiles_img(smiles, border="1"))
            if len(row) == ROW_SIZE:
                imgs.append(ui.div(row))
                row = []
        if row:
            imgs.append(ui.div(row))
        return ui.div(imgs)


www_dir = Path(__file__).parent / "www"
www_dir.mkdir(parents=True, exist_ok=True)
app = App(app_ui, server, static_assets=www_dir)
