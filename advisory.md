# PA-app4 Advisory

## Known issues and recommendations

### 1. uv must be installed separately
`pyproject.toml` assumes `uv` is available. Install with:
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```
Then from `backend/`:
```bash
uv sync --extra dev
```
If `uv` is unavailable, fall back to:
```bash
python -m venv venv && source venv/bin/activate && pip install -e ".[dev]"
```

### 2. Astro requires Node ≥ 18
Ensure Node.js ≥ 18 is installed. Check with `node --version`.

### 3. SVG imports in React (Astro)
`gecko-investigate.svg` is imported in `GeckoInvestigate.tsx` as a static asset. Astro returns an object `{ src: string }` for image imports rather than a plain string. The component handles this with `geckoIcon.src ?? geckoIcon`. If you see a broken image, verify Astro's image handling for SVGs in your version.

### 4. Prot_modules.py / intron_modules.py — no type stubs
These large modules are not typed and are excluded from strict mypy checking via `ignore_missing_imports = true` in `pyproject.toml`. Adding type stubs would require significant effort and is not recommended until the modules stabilise.

### 5. Tests require the full Python environment
`test_meta.py` imports `feature_detection.py` which imports `biopython`. The full `uv sync` must be run before tests will pass. Run tests with:
```bash
cd backend
uv run pytest tests/ -q
```

### 6. No PostgreSQL / Supabase integration yet
CLAUDE.md recommends Postgres/Supabase for data storage. The current app stores nothing persistently — all processing is stateless. This is appropriate for the current use case. If result persistence is needed in future, consider adding a Supabase client.
