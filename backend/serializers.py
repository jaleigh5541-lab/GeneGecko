from dataclasses import asdict
from typing import Any


def serialize_results(results: list[Any]) -> list[dict[str, Any]]:
    """Convert processing results to JSON-serializable dicts."""
    out = []
    for r in results:
        r = dict(r)
        if r.get("detected_features"):
            r["detected_features"] = [asdict(f) for f in r["detected_features"]]
        out.append(r)
    return out
