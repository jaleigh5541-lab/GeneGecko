from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from routes.process import router as process_router
from routes.meta import router as meta_router
from routes.align import router as align_router

app = FastAPI(title="Gene-Gecko API")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(process_router)
app.include_router(meta_router)
app.include_router(align_router)

if __name__ == "__main__":
    import uvicorn
    uvicorn.run("main:app", host="0.0.0.0", port=5000, reload=True)
