from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Iterable


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def iter_files(root: Path) -> Iterable[Path]:
    for path in sorted(root.rglob("*")):
        if path.is_file() and ".git" not in path.parts:
            yield path


def write_manifest(root: str | Path, out_dir: str | Path) -> Path:
    root = Path(root).resolve()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    manifest = out / "files_manifest.csv"
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["path", "size_bytes", "sha256"])
        writer.writeheader()
        for path in iter_files(root):
            writer.writerow(
                {
                    "path": str(path.relative_to(root)).replace("\\", "/"),
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
            )
    return manifest


def write_ro_crate(root: str | Path, out_path: str | Path) -> Path:
    root = Path(root).resolve()
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    graph = [
        {
            "@id": "ro-crate-metadata.json",
            "@type": "CreativeWork",
            "about": {"@id": "./"},
        },
        {
            "@id": "./",
            "@type": "Dataset",
            "name": "ORBIT benchmark run",
            "description": "Files produced or consumed by an ORBIT benchmark run.",
            "hasPart": [{"@id": str(p.relative_to(root)).replace("\\", "/")} for p in iter_files(root)],
        },
    ]
    payload = {"@context": "https://w3id.org/ro/crate/1.1/context", "@graph": graph}
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return out

