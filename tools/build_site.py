#!/usr/bin/env python3
"""Build the research homepage, MkDocs documentation and source download."""
from __future__ import annotations

import gzip
import json
from pathlib import Path
import shutil
import zipfile
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "website/dist"
DIRECTORIES = ("src", "examples", "integrations", "data", "test", "analysis", "demo", "docs", "tools", "website", ".github")
FILES = ("README.md", "Project.toml", "CITATION.cff", ".gitignore", "mkdocs.yml", "requirements-docs.txt", "requirements-figures.txt", "references/README.md",
         "docs/assets/publications/FSchnell_BAThesis_AMMicrostructure.pdf")


def source_files():
    """Use explicit roots; never traverse local references, environments or output."""
    paths = [ROOT / name for name in FILES]
    for name in DIRECTORIES:
        for path in (ROOT / name).rglob("*"):
            rel = path.relative_to(ROOT)
            if path.is_file() and not path.is_symlink() and not any(
                part in {"dist", "__pycache__", "node_modules", ".DS_Store"} for part in rel.parts
            ) and path.suffix.lower() not in {".pdf", ".deb", ".pyc"}:
                paths.append(path)
    return sorted({p for p in paths if p.is_file()}, key=lambda p: p.relative_to(ROOT).as_posix())


def build():
    files = source_files()  # Capture before creating the generated directory.
    if DEST.exists():
        shutil.rmtree(DEST)
    DEST.mkdir(parents=True)
    for name in ("index.html", "styles.css", "app.js"):
        shutil.copy2(ROOT / "website" / name, DEST / name)
    shutil.copytree(ROOT / "website/assets", DEST / "assets")
    shutil.copytree(ROOT / "docs/assets", DEST / "assets", dirs_exist_ok=True)
    subprocess.run([sys.executable, "-m", "mkdocs", "build", "--strict", "--config-file", str(ROOT / "mkdocs.yml")], check=True)
    data = json.loads(gzip.decompress((ROOT / "docs/assets/demo/grains.json.gz").read_bytes()))
    (DEST / "assets/demo-data.js").write_text("window.AM_DEMO=" + json.dumps(data, separators=(",", ":"), allow_nan=False) + ";\n")
    timeline = json.loads((ROOT / "docs/assets/demo/evolution.json").read_text())
    (DEST / "assets/evolution-data.js").write_text("window.AM_EVOLUTION=" + json.dumps(timeline, separators=(",", ":"), allow_nan=False) + ";\n")
    (DEST / "downloads").mkdir()
    with zipfile.ZipFile(DEST / "downloads/ammicrostructure-source.zip", "w", zipfile.ZIP_DEFLATED) as archive:
        for path in files:
            archive.write(path, "AM_Microstructure/" + path.relative_to(ROOT).as_posix())
    (DEST / ".nojekyll").touch()
    print(f"Built {DEST} ({len(files)} source files in download)")


if __name__ == "__main__":
    build()
