#!/bin/bash
set -e

ROOT="$(cd "$(dirname "$0")" && pwd)"
BACKEND="$ROOT/backend"
FRONTEND="$ROOT/frontend"

# Start FastAPI backend with uvicorn (using uv if available, else direct uvicorn)
echo "Starting FastAPI backend on http://localhost:5000 ..."
cd "$BACKEND"
source "$HOME/pa-app4-venv/bin/activate"
python main.py &
BACKEND_PID=$!

# Start Astro frontend
echo "Starting Astro frontend ..."
cd "$FRONTEND"
npm run dev &
FRONTEND_PID=$!

echo ""
echo "App running:"
echo "  Backend:  http://localhost:5000"
echo "  Frontend: http://localhost:5173"
echo ""
echo "Press Ctrl+C to stop both servers."

trap "echo 'Stopping...'; kill $BACKEND_PID $FRONTEND_PID 2>/dev/null; exit 0" SIGINT SIGTERM

wait
