import gc
import json
import re
import threading
import time
import uuid
from typing import Any, Optional

from fastapi import APIRouter, File, Form, HTTPException, UploadFile
from fastapi.responses import JSONResponse

from adapters import FileWrapper
from serializers import serialize_results

router = APIRouter(prefix="/api")

# In-memory job store
_jobs: dict[str, dict[str, Any]] = {}
_JOBS_TTL = 3600


def _cleanup_jobs() -> None:
    now = time.time()
    expired = [jid for jid, j in _jobs.items() if now - j["created_at"] > _JOBS_TTL]
    for jid in expired:
        del _jobs[jid]


def _process_paired_one_at_a_time(
    job_id: str,
    seq_wrappers: list[FileWrapper],
    ref_wrapper: FileWrapper,
    use_dna_alignment: bool,
    min_lcs: int,
    selected: Optional[list[str]],
) -> list[dict[str, Any]]:
    """Process T7/T7-Term pairs one at a time to minimize memory."""
    from Prot_modules import (  # type: ignore[import]
        find_pairs,
        process_pair,
        load_reference_xlsx,
        load_reference_dna_xlsx,
        detect_features,
        _add_orf_info,
    )

    filenames = [w.name for w in seq_wrappers]
    name_to_wrapper = {w.name: w for w in seq_wrappers}

    pairs = find_pairs(filenames)
    if not pairs:
        raise ValueError(
            "No valid T7/T7-Term pairs found. "
            "Files should follow naming pattern: *-T7.seq and *-T7-Term.seq"
        )

    # Load references once
    if use_dna_alignment:
        ref_dna_dict = load_reference_dna_xlsx(ref_wrapper)
        ref_dict = None
    else:
        ref_dict = load_reference_xlsx(ref_wrapper)
        ref_dna_dict = None

    total_pairs = len(pairs)
    all_results: list[dict[str, Any]] = []

    for idx, (fwd_name, rev_name, prefix) in enumerate(pairs):
        fwd_file = name_to_wrapper[fwd_name]
        rev_file = name_to_wrapper[rev_name]

        try:
            # Reset file positions
            fwd_file.seek(0)
            rev_file.seek(0)
            ref_wrapper.seek(0)

            result = process_pair(
                fwd_file, rev_file, prefix, min_lcs=min_lcs,
                ref_dict=ref_dict, ref_dna_dict=ref_dna_dict,
            )
            result["sample"] = prefix
            result["pipeline"] = "protein"
            result["orientation"] = "forward"
            result["detected_features"] = detect_features(
                result.get("merged_dna") or "", selected_categories=selected
            ) if result.get("merged_dna") else []
            result["crop_info"] = {}
            result["final_dna"] = result.get("merged_dna")

            sample_base = re.split(r'[_-]', prefix, 1)[0].lower()
            if ref_dna_dict:
                result["reference"] = ref_dna_dict.get(sample_base)
            elif ref_dict:
                result["reference"] = ref_dict.get(sample_base)

            _add_orf_info(result, result.get("merged_dna"))

        except Exception as e:
            result = {
                "sample": prefix, "error": str(e),
                "merged_dna": None, "final_dna": None,
                "longest_orf": None, "orf_dna": None, "orf_length": 0,
                "frame": None, "has_orf": False,
                "merge_info": {}, "crop_info": {},
                "reference": None, "alignment_text": None,
                "identity_percent": 0.0, "orientation": "forward",
                "pipeline": "protein", "detected_features": [],
            }

        serialized = serialize_results([result])
        all_results.extend(serialized)
        del result, serialized
        gc.collect()

        _jobs[job_id]["progress"] = {
            "processed": idx + 1,
            "total": total_pairs,
        }

    return all_results


def _process_batch(
    seq_batch: list[FileWrapper],
    ref_wrapper: FileWrapper,
    sequence_type: str,
    use_protein: bool,
    use_dna_alignment: bool,
    min_lcs: int,
    selected: Optional[list[str]],
) -> list[dict[str, Any]]:
    """Process a batch for non-paired modes."""
    ref_wrapper.seek(0)
    if use_protein:
        from Prot_modules import (  # type: ignore[import]
            process_circular_queries as prot_circular,
            process_unidirectional as prot_unidirectional,
        )
        if sequence_type == "circular":
            results = prot_circular(
                seq_batch, ref_wrapper,
                use_dna_alignment=use_dna_alignment,
                selected_categories=selected,
            )
        else:
            results = prot_unidirectional(
                seq_batch, ref_wrapper,
                use_dna_alignment=use_dna_alignment,
                selected_categories=selected,
            )
    else:
        from intron_modules import (  # type: ignore[import]
            process_circular_queries as intron_circular,
            process_unidirectional as intron_unidirectional,
        )
        if sequence_type == "circular":
            results = intron_circular(
                seq_batch, ref_wrapper,
                selected_categories=selected,
            )
        else:
            results = intron_unidirectional(
                seq_batch, ref_wrapper,
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

        if sequence_type == "paired":
            # Process one pair at a time — minimal memory
            all_results = _process_paired_one_at_a_time(
                job_id, seq_wrappers, ref_wrapper,
                use_dna_alignment, min_lcs, selected,
            )
        else:
            # Non-paired: batch 3 files at a time
            batch_size = 3
            total = len(seq_wrappers)
            all_results = []

            for i in range(0, total, batch_size):
                batch = seq_wrappers[i : i + batch_size]
                for w in batch:
                    w.seek(0)
                batch_results = _process_batch(
                    batch, ref_wrapper, sequence_type,
                    use_protein, use_dna_alignment, min_lcs, selected,
                )
                all_results.extend(batch_results)

                _jobs[job_id]["progress"] = {
                    "processed": min(i + batch_size, total),
                    "total": total,
                }

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
