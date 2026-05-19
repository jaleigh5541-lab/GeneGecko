"""Tests for /api/health and /api/categories endpoints.

These tests are written first (TDD) and define the expected API contract.
"""
import sys
import os

# Allow imports from the backend root
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import pytest
from httpx import AsyncClient, ASGITransport

from main import app


@pytest.mark.asyncio
async def test_health_returns_ok() -> None:
    async with AsyncClient(
        transport=ASGITransport(app=app), base_url="http://test"
    ) as client:
        response = await client.get("/api/health")
    assert response.status_code == 200
    assert response.json() == {"status": "ok"}


@pytest.mark.asyncio
async def test_categories_returns_expected_shape() -> None:
    async with AsyncClient(
        transport=ASGITransport(app=app), base_url="http://test"
    ) as client:
        response = await client.get("/api/categories")
    assert response.status_code == 200
    data = response.json()
    assert "categories" in data
    assert "colors" in data
    assert "orf_min_aa" in data
    assert isinstance(data["categories"], list)
    assert isinstance(data["colors"], dict)
