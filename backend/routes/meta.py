from fastapi import APIRouter
from fastapi.responses import JSONResponse

from feature_detection import get_all_categories, CATEGORY_COLORS  # type: ignore[import]

router = APIRouter(prefix="/api")


@router.get("/categories")
def categories() -> JSONResponse:
    return JSONResponse({
        "categories": get_all_categories(),
        "colors": CATEGORY_COLORS,
        "orf_min_aa": 30,
    })


@router.get("/health")
def health() -> JSONResponse:
    return JSONResponse({"status": "ok"})
