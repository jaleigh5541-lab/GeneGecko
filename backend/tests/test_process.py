"""Tests for /api/process endpoint.

These tests are written first (TDD) and define the expected API contract.
Integration tests requiring full processing modules are skipped in CI.
"""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import pytest
from httpx import AsyncClient, ASGITransport

from main import app


@pytest.mark.asyncio
async def test_process_missing_files_returns_422() -> None:
    """POST /api/process with no files should return 422 (validation error)."""
    async with AsyncClient(
        transport=ASGITransport(app=app), base_url="http://test"
    ) as client:
        response = await client.post("/api/process", data={})
    assert response.status_code == 422


@pytest.mark.asyncio
async def test_process_with_dummy_files_returns_json() -> None:
    """POST /api/process with minimal dummy files should return JSON (may error
    in processing, but the API layer should respond with JSON)."""
    dummy_seq = b">sample1\nATGCATGC\n"
    dummy_ref = b"dummy xlsx content"

    async with AsyncClient(
        transport=ASGITransport(app=app), base_url="http://test"
    ) as client:
        response = await client.post(
            "/api/process",
            files={
                "seq_files": ("sample1_T7.seq", dummy_seq, "application/octet-stream"),
                "ref_file": ("reference.xlsx", dummy_ref, "application/octet-stream"),
            },
            data={"sequence_type": "paired", "primer_type": "protein"},
        )
    # Response should always be JSON (success or error detail)
    assert response.headers["content-type"].startswith("application/json")
