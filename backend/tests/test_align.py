"""Tests for /api/align endpoint."""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import pytest
from httpx import AsyncClient, ASGITransport

from main import app


@pytest.mark.asyncio
async def test_align_dna_returns_alignment() -> None:
    async with AsyncClient(
        transport=ASGITransport(app=app), base_url="http://test"
    ) as client:
        response = await client.post(
            "/api/align",
            json={
                "seq1": "ATGATGATGATGATGATG",
                "seq2": "ATGATGATGATGATGATG",
                "seq_type": "dna",
            },
        )
    assert response.status_code == 200
    data = response.json()
    assert data["success"] is True
    assert "alignment_html" in data
    assert data["identity_percent"] == 100.0


@pytest.mark.asyncio
async def test_align_protein_returns_alignment() -> None:
    async with AsyncClient(
        transport=ASGITransport(app=app), base_url="http://test"
    ) as client:
        response = await client.post(
            "/api/align",
            json={
                "seq1": "MFILTERPROTEIN",
                "seq2": "MFILTERPROTEIN",
                "seq_type": "protein",
            },
        )
    assert response.status_code == 200
    data = response.json()
    assert data["success"] is True
    assert data["identity_percent"] == 100.0
    assert data["longest_orf_length"] > 0


@pytest.mark.asyncio
async def test_align_short_sequences_rejected() -> None:
    async with AsyncClient(
        transport=ASGITransport(app=app), base_url="http://test"
    ) as client:
        response = await client.post(
            "/api/align",
            json={"seq1": "AT", "seq2": "AT", "seq_type": "dna"},
        )
    assert response.status_code == 400


@pytest.mark.asyncio
async def test_align_dna_orf_detection() -> None:
    # A sequence containing an ATG start codon should produce ORFs
    seq = "ATGAAAGGGCCCTTTTAG"  # M K G P F *
    async with AsyncClient(
        transport=ASGITransport(app=app), base_url="http://test"
    ) as client:
        response = await client.post(
            "/api/align",
            json={"seq1": seq, "seq2": seq, "seq_type": "dna"},
        )
    assert response.status_code == 200
    data = response.json()
    assert data["orf_count"] >= 1
    assert data["longest_orf_length"] >= 1
