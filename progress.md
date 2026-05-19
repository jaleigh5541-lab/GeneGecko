# PA-app4 Progress

## 2026-03-17 — Initial rebuild per CLAUDE.md

### Completed

**Backend (Flask → FastAPI)**
- Replaced `app.py` (Flask) with `main.py` (FastAPI + uvicorn)
- Rewrote `routes/process.py` using `APIRouter`, `Form(...)`, `UploadFile`, async file reading
- Rewrote `routes/meta.py` using `APIRouter`
- Updated `adapters.py`: `FileWrapper` now takes `(filename: str, data: bytes)` directly instead of a Flask `FileStorage` object
- Updated `serializers.py` with static type annotations
- Created `pyproject.toml` for uv dependency management (replaces `requirements.txt`)
  - Dependencies: fastapi, uvicorn[standard], python-multipart, biopython, pandas, openpyxl
  - Dev deps: pytest, httpx, mypy, ruff, pytest-asyncio
- Kept compute modules as-is: `Prot_modules.py`, `intron_modules.py`, `feature_detection.py`

**Tests (TDD)**
- `tests/test_meta.py`: async tests for `/api/health` and `/api/categories`
- `tests/test_process.py`: validation tests for `/api/process` (missing files → 422, dummy files → JSON response)

**Frontend (Vite+React JSX → Astro + React + Tailwind + TypeScript)**
- New Astro project in `frontend/`
- `astro.config.mjs` with `@astrojs/react` and `@astrojs/tailwind` integrations
- API proxy `/api → http://localhost:5000` via Vite config
- `tailwind.config.mjs`, `tsconfig.json`, `.eslintrc.cjs`
- All JSX components ported to TSX with typed props
- `src/api/client.ts` — typed API client
- `src/hooks/useLoadingText.ts` — typed loading text hook
- `src/styles/global.css` — Tailwind directives + full Gene-Gecko custom CSS
- `src/pages/index.astro` — Astro page that mounts `<App client:load />`
- Assets copied: `lilthang.json`, `toungy-lizard.json`, `gecko-investigate.svg`, `favicon.svg`

**Run script**
- `run.sh` updated to use `uvicorn` (via `uv run` if available) and `npm run dev` (Astro)

### To do before first run
1. `cd backend && uv sync --extra dev` — install Python deps
2. `cd frontend && npm install` — install Node deps
3. `./run.sh` — start both servers
