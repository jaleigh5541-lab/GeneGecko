import json
from typing import Optional

from fastapi import APIRouter, File, Form, HTTPException, UploadFile
from fastapi.responses import JSONResponse

from adapters import FileWrapper
from serializers import serialize_results

router = APIRouter(prefix="/api")


@router.post("/process")
async def process(
    seq_files: list[UploadFile] = File(...),
    ref_file: UploadFile = File(...),
    sequence_type: str = Form("paired"),
    primer_type: str = Form("auto"),
    use_aa_refs: str = Form("false"),
    min_lcs: int = Form(12),
    selected_categories: Optional[str] = Form(None),
) -> JSONResponse:
    if not seq_files:
        raise HTTPException(status_code=400, detail="No sequence files uploaded.")

    selected: Optional[list[str]] = (
        json.loads(selected_categories) if selected_categories else None
    )

    seq_wrappers = [
        FileWrapper(f.filename or "", await f.read()) for f in seq_files
    ]
    ref_wrapper = FileWrapper(ref_file.filename or "", await ref_file.read())

    use_protein = primer_type != "intron"
    use_dna_alignment = use_aa_refs.lower() != "true"

    try:
        if use_protein:
            from Prot_modules import (  # type: ignore[import]
                main as prot_main,
                process_circular_queries as prot_circular,
                process_unidirectional as prot_unidirectional,
            )

            if sequence_type == "circular":
                results = prot_circular(
                    seq_wrappers, ref_wrapper,
                    use_dna_alignment=use_dna_alignment,
                    selected_categories=selected,
                )
            elif sequence_type == "unidirectional":
                results = prot_unidirectional(
                    seq_wrappers, ref_wrapper,
                    use_dna_alignment=use_dna_alignment,
                    selected_categories=selected,
                )
            else:
                results = prot_main(
                    seq_wrappers, ref_wrapper,
                    min_lcs=min_lcs,
                    use_dna_alignment=use_dna_alignment,
                )
        else:
            from intron_modules import (  # type: ignore[import]
                main as intron_main,
                process_circular_queries as intron_circular,
                process_unidirectional as intron_unidirectional,
            )

            if sequence_type == "circular":
                results = intron_circular(
                    seq_wrappers, ref_wrapper,
                    selected_categories=selected,
                )
            elif sequence_type == "unidirectional":
                results = intron_unidirectional(
                    seq_wrappers, ref_wrapper,
                    selected_categories=selected,
                )
            else:
                results = intron_main(
                    seq_wrappers, ref_wrapper,
                    min_lcs=min_lcs,
                    selected_categories=selected,
                )

        return JSONResponse({"success": True, "results": serialize_results(results)})

    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e)) from e
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e)) from e
