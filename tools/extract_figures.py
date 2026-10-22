#!/usr/bin/env python3
"""Recreate the documented research figure crops from the local source PDFs.

Requires Python 3 and Poppler's ``pdftoppm``; no Python packages are needed.
The thesis is included in docs/assets/publications/. Place the journal article
in references/local/ before running this script.

Usage: python3 tools/extract_figures.py [--dpi 240] [--only paper-pipeline]
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import struct
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs" / "assets" / "research"

SOURCES = {
    "paper": {
        "filename": "1-s2.0-S026412752601289X-main.pdf",
        "title": "Experimentally informed physics-based multi-scale simulation for microstructure tailoring and cracking behavior modification in Inconel 718 during PBF-LB/M",
        "authors": "Hannes Panzer, Saeed Maleksaeedi, Hamidreza Aghajani, Mohsen K. Keshavarz, Fangzhou Li, Ferdinand Schnell, Kai-Uwe Beuerlein, Felix Riegger, Michael F. Zaeh",
        "year": 2026,
        "publication": "Materials & Design 269, 116714",
        "url": "https://doi.org/10.1016/j.matdes.2026.116714",
        "license": "CC BY 4.0",
        "license_url": "https://creativecommons.org/licenses/by/4.0/",
        "license_evidence": "Source PDF page 1: copyright 2026 The Authors; published by Elsevier Ltd.; open access under CC BY 4.0.",
        "attribution": "Panzer et al. (2026), Materials & Design 269, 116714, doi:10.1016/j.matdes.2026.116714. CC BY 4.0. Cropped from the published figure; content unchanged.",
    },
    "thesis": {
        "filename": "FSchnell_BAThesis_AMMicrostructure.pdf",
        "title": "Simulating Grain Growth in PBF-LB/M with a Cellular Automaton Model to Determine Fracture Criticality",
        "authors": "Ferdinand O. Schnell",
        "year": 2024,
        "publication": "Bachelor's thesis, Technical University of Munich, submitted 31 May 2024",
        "url": None,
        "license": "Copyright Ferdinand O. Schnell; thesis and figure excerpts included with the author's permission. No general reuse license declared.",
        "license_url": None,
        "license_evidence": "The author requested inclusion of the full thesis and figure excerpts in this repository and website. The PDF does not declare a Creative Commons license.",
        "attribution": "Ferdinand O. Schnell (2024), bachelor's thesis, Technical University of Munich. Cropped from the thesis; content unchanged. No general reuse license declared.",
    },
}

# Crops are PDF points (1/72 inch), measured from the page's top-left corner.
# Each crop includes the complete chosen figure, its axes, legends, and panel
# labels, while omitting page furniture and the original typeset caption. Captions
# and provenance remain available in the generated manifest and research.md.
FIGURES = [
    {
        "id": "paper-pipeline",
        "source": "paper",
        "figure": "2",
        "pdf_page": 4,
        "printed_page": "4",
        "crop_points": {"x": 116, "y": 56, "width": 364, "height": 220},
        "caption": "The published multi-scale workflow: a global crack analysis identifies critical regions, refined thermal simulations drive local grain growth, and microstructural descriptors feed back into the failure model.",
        "original_caption": "Approach for the numerical prediction of cracking as a function of the microstructure; z is the building direction.",
        "alt": "Four-stage workflow from macro-scale crack-critical regions through meso-scale thermal fields and micro-scale grain simulation to a microstructure-dependent cracking model.",
        "scope_note": "Wider research pipeline in the 2026 paper, not a claim that every stage is provided by this repository or solely authored by Ferdinand Schnell.",
    },
    {
        "id": "paper-microstructure-comparison",
        "source": "paper",
        "figure": "10",
        "pdf_page": 10,
        "printed_page": "10",
        "crop_points": {"x": 87, "y": 64, "width": 416, "height": 416},
        "caption": "As-built microstructures in the later study: measured EBSD maps and pole figures for PP0 and PP2 above, simulated grain maps below. Grain morphology improved, while the predicted crystallographic texture remained imperfect.",
        "original_caption": "Microstructure characterization in the as-built state: (a) EBSD image (PP0), (b) pole figures (PP0), (c) EBSD image (PP2), (d) pole figures (PP2), (e) simulated EBSD image (PP0), (f) simulated EBSD image (PP2); z is the building direction.",
        "alt": "Measured colored EBSD maps and circular pole figures for PP0 and PP2, with two simulated colored grain maps below; scale bars are 100 micrometres.",
        "scope_note": "The 2026 paper's model and results; not a fresh run of the thesis code or the repository's small example.",
    },
    {
        "id": "paper-grain-statistics",
        "source": "paper",
        "figure": "11",
        "pdf_page": 11,
        "printed_page": "11",
        "crop_points": {"x": 115, "y": 55, "width": 368, "height": 148},
        "caption": "Measured and simulated grain characteristics across PP0-PP3. The later model captured increasing grain size and aspect ratio, with reported maximum deviations of 20% and 16%; it did not reproduce the experimental Schmid-factor trend.",
        "original_caption": "Comparison of the experimentally and simulatively determined grain characteristics: (a) grain sizes and aspect ratios, (b) Schmid factors.",
        "alt": "Two charts compare experiment and simulation across process parameter sets PP0 to PP3: grain size and aspect ratio on the left, Schmid factor on the right.",
        "scope_note": "The 16% aspect-ratio deviation here is a different quantity from the paper's approximately 16% reduction in mean crack length.",
    },
    {
        "id": "paper-crack-lengths",
        "source": "paper",
        "figure": "13",
        "pdf_page": 12,
        "printed_page": "12",
        "crop_points": {"x": 45, "y": 55, "width": 237, "height": 185},
        "caption": "Experimentally measured crack lengths for four process parameter sets. The paper reports an approximately 16% shorter mean crack length in the representative PP3-versus-PP0 crack-shape comparison discussed with Fig. 15 on the same page; this is a collaborative experimental result.",
        "original_caption": "Experimentally identified crack lengths with changing process parameters.",
        "alt": "Box plots of mean crack length in millimetres for PP0, PP1, PP2, and PP3; red crosses and a dashed line mark the means.",
        "scope_note": "Experimental process-parameter comparison in the 2026 paper. The 16% statement is reported in the discussion of Fig. 15(b), PP3 (350 W, no hatch rotation), compared with PP0 (285 W, rotating hatch); it is not a percentage calculated from this crop, a normalized simulated crack length, or a performance claim for this repository's standalone example.",
    },
    {
        "id": "thesis-single-track",
        "source": "thesis",
        "figure": "4.1",
        "pdf_page": 56,
        "printed_page": "42",
        "crop_points": {"x": 102, "y": 418, "width": 410, "height": 90},
        "caption": "A 2024 thesis verification case: with bulk nucleation disabled (N0 = 0), a single melt track solidifies by growth from existing grains at the melt-pool boundary.",
        "original_caption": "Single melt track with nucleation N0 = 0, where all growth stems from existing grains along the melt pool boundary.",
        "alt": "A long, narrow, colored single melt-track cross section shows grains growing from an existing melt-pool boundary; a color scale shows grain orientation angle.",
        "scope_note": "A controlled verification result, not evidence that the thesis model reproduced realistic multi-layer microstructure.",
    },
    {
        "id": "thesis-multilayer",
        "source": "thesis",
        "figure": "4.7",
        "pdf_page": 60,
        "printed_page": "46",
        "crop_points": {"x": 90, "y": 116, "width": 420, "height": 240},
        "caption": "A revealing thesis result: grains can cross the yellow build-layer boundaries, but irregular nucleation prevents realistic columnar growth and consistent grain-size trends. The circles mark contrasting local grain structures in a 285 W case.",
        "original_caption": "Grain structure in the X-Z plane cut in the middle, with the powder bed height at each layer marked in yellow. Output of Simulation Scenario 2 with 285 W laser power.",
        "alt": "Colored grain structure in an X-Z cross section with horizontal yellow layer boundaries and three circles highlighting differently sized and oriented grains.",
        "scope_note": "Read with thesis printed pages 45-47 and 49-51: this is evidence of both implemented multi-layer growth and the model's validation limitations.",
    },
]


def locate_source(filename: str) -> Path:
    for folder in (ROOT / "docs/assets/publications", ROOT / "references" / "local", ROOT / "references", ROOT):
        candidate = folder / filename
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(
        f"Missing {filename}. Place the locally held source PDF in references/local/."
    )


def png_size(path: Path) -> tuple[int, int]:
    with path.open("rb") as stream:
        header = stream.read(24)
    if header[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError(f"Invalid PNG output: {path}")
    return struct.unpack(">II", header[16:24])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dpi", type=int, default=240)
    parser.add_argument("--only", choices=[figure["id"] for figure in FIGURES])
    args = parser.parse_args()
    if not 72 <= args.dpi <= 600:
        parser.error("--dpi must be between 72 and 600")
    if not shutil.which("pdftoppm"):
        parser.error("pdftoppm is required (install the Poppler utilities)")
    OUTPUT.mkdir(parents=True, exist_ok=True)

    selected = [f for f in FIGURES if not args.only or f["id"] == args.only]
    source_paths = {key: locate_source(SOURCES[key]["filename"]) for key in {f["source"] for f in selected}}
    generated = []
    for figure in selected:
        crop = figure["crop_points"]
        pixels = {key: round(value * args.dpi / 72) for key, value in crop.items()}
        prefix = OUTPUT / figure["id"]
        subprocess.run(
            [
                "pdftoppm", "-f", str(figure["pdf_page"]), "-l", str(figure["pdf_page"]),
                "-singlefile", "-r", str(args.dpi), "-x", str(pixels["x"]),
                "-y", str(pixels["y"]), "-W", str(pixels["width"]),
                "-H", str(pixels["height"]), "-png", str(source_paths[figure["source"]]),
                str(prefix),
            ],
            check=True,
        )
        path = prefix.with_suffix(".png")
        width, height = png_size(path)
        generated.append({
            **figure,
            "file": path.name,
            "dpi": args.dpi,
            "crop_pixels": pixels,
            "width": width,
            "height": height,
            "bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "transformations": "Rasterized with Poppler and cropped; original figure content, colors, labels, and annotations retained. Original typeset caption omitted from image and recorded here.",
        })
        print(f"{path.relative_to(ROOT)}: {width} x {height}")

    manifest_path = OUTPUT / "manifest.json"
    if args.only and manifest_path.exists():
        previous = json.loads(manifest_path.read_text())["figures"]
        replacements = {figure["id"]: figure for figure in generated}
        generated = [replacements.get(figure["id"], figure) for figure in previous]
        existing = {figure["id"] for figure in generated}
        generated.extend(figure for key, figure in replacements.items() if key not in existing)
    sources = {}
    for key, source in SOURCES.items():
        record = dict(source)
        if key in source_paths:
            record["sha256"] = hashlib.sha256(source_paths[key].read_bytes()).hexdigest()
        elif manifest_path.exists():
            record.update(json.loads(manifest_path.read_text())["sources"].get(key, {}))
        sources[key] = record
    manifest = {
        "schema_version": 1,
        "generator": "tools/extract_figures.py",
        "coordinate_system": "PDF points (1/72 inch) from top-left; PDF page numbers are 1-based. Printed page labels are separate.",
        "sources": sources,
        "figures": generated,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, subprocess.CalledProcessError) as error:
        print(error, file=sys.stderr)
        raise SystemExit(1)
