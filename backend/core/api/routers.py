import io
import asyncio
import logging

from fastapi import APIRouter
from fastapi.responses import StreamingResponse

from api.services.draw import draw_molecule

log = logging.getLogger(__name__)

router = APIRouter()


@router.post("/substructure/")
@router.post("/substructure", include_in_schema=False)
async def substructure(
        smiles: str
):
    await asyncio.sleep(1)
    result = [smiles]
    return result


@router.get("/render/")
@router.get("/render", include_in_schema=False)
async def render(
        q: str,
        height: int = 450,
        width: int = 150,
        r: float = 1,
        g: float = 1,
        b: float = 1,
):
    res: io.StringIO = draw_molecule(q, r, g, b, molSize=(height, width))

    block_size = 4096
    value = iter(lambda: res.read(block_size), "")
    return StreamingResponse(value, media_type="image/svg+xml")
