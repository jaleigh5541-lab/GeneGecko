import gc
import json
import threading
import time
import uuid
from typing import Any, Optional

from fastapi import APIRouter, File, Form, HTTPException, UploadFile
from fastapi.responses import JSONResponse

from adapters import FileWrapper
from serializers import serialize_results

router = APIRouter(prefix="/api")

BATCH_SIZE = 5  # process this many files at a time to limit peak memory

# In-memory job store: {job_id: {status, results, error, created_at, progress}}
_jobs: dict[str, dict[str, Any]] = {}
_JOBS_TTL = 3600


def _cleanup_jobs() -> None:
    now = time.time()
    expired = [jid for jid, j in _jobs.items() if now - j["created_at"] > _JOBS_TTL]
    for jid in expired:
        del _jobs[jid]


def _process_batch(
    seq_batch: list[FileWrapper],
    ref_wrapper: FileWrapper,
    sequence_type: str,
    use_protein: bool,
    use_dna_alignment: bool,
    min_lcs: int,
    selected: Optional[list[str]],
) -> list[dict[str, Any]]:
    """Process a single batch of sequence files and return serialized results."""
    if use_protein:
        from Prot_modules import (  # type: ignore[import]
            main as prot_main,
            process_circular_queries as prot_circular,
            process_unidirectional as prot_unidirectional,
        )
        if sequence_type == "circular":
            results = prot_circular(
                seq_batch, ref_wrapper,
                use_dna_alignment=use_dna_alignment,
                selected_categories=selected,
            )
        elif sequence_type == "unidirectional":
            results = prot_unidirectional(
                seq_batch, ref_wrapper,
                use_dna_alignment=use_dna_alignment,
                selected_categories=selected,
            )
        else:
            results = prot_main(
                seq_batch, ref_wrapper,
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
                seq_batch, ref_wrapper,
                selected_categories=selected,
            )
        elif sequence_type == "unidirectional":
            results = intron_unidirectional(
                seq_batch, ref_wrapper,
                selected_categories=selected,
            )
        else:
            results = intron_main(
                seq_batch, ref_wrapper,
                min_lcs=min_lcs,
                selected_categories=selected,
            )

    serialized = serialize_results(results)
    del results
    gc.collect()
    return serialized


def _run_processing(
    job_id: str,
    seq_wrappers: list[FileWrapper],
    ref_wrapper: FileWrapper,
    sequence_type: str,
    primer_type: str,
    use_aa_refs: str,
    min_lcs: int,
    selected: Optional[list[str]],
) -> None:
    try:
        use_protein = primer_type != "intron"
        use_dna_alignment = use_aa_refs.lower() != "true"
        total = len(seq_wrappers)

        all_results: list[dict[str, Any]] = []

        for i in range(0, total, BATCH_SIZE):
            batch = seq_wrappers[i : i + BATCH_SIZE]
            batch_results = _process_batch(
                batch, ref_wrapper, sequence_type,
                use_protein, use_dna_alignment, min_lcs, selected,
            )
            all_results.extend(batch_results)

            # Update progress
            processed = min(i + BATCH_SIZE, total)
            _jobs[job_id]["progress"] = {
                "processed": processed,
                "total": total,
            }

            # Free the batch
            del batch, batch_results
            gc.collect()

        _jobs[job_id]["results"] = all_results
        _jobs[job_id]["status"] = "done"

    except Exception as e:
        _jobs[job_id]["status"] = "error"
        _jobs[job_id]["error"] = str(e)


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

    _cleanup_jobs()

    job_id = uuid.uuid4().hex
    _jobs[job_id] = {
        "status": "running",
        "results": None,
        "error": None,
        "created_at": time.time(),
        "progress": {"processed": 0, "total": len(seq_wrappers)},
    }

    thread = threading.Thread(
        target=_run_processing,
        args=(job_id, seq_wrappers, ref_wrapper, sequence_type, primer_type,
              use_aa_refs, min_lcs, selected),
        daemon=True,
    )
    thread.start()

    return JSONResponse({"success": True, "job_id": job_id})


@router.get("/jobs/{job_id}")
async def get_job(job_id: str) -> JSONResponse:
    job = _jobs.get(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    if job["status"] == "running":
        return JSONResponse({
            "status": "running",
            "progress": job.get("progress"),
        })
    elif job["status"] == "error":
        return JSONResponse({"status": "error", "error": job["error"]})
    else:
        return JSONResponse({
            "status": "done",
            "success": True,
            "results": job["results"],
        })
