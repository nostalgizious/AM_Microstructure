#!/usr/bin/env python3
"""Export real VTK cell fields to the checked-in static viewer data; stdlib only."""
import argparse
import csv
import colorsys
import gzip
import hashlib
import json
import math
from pathlib import Path
import tomllib

ROOT = Path(__file__).resolve().parents[1]


def write_gzip_json(path, data):
    """Deterministic compression keeps full-resolution scientific fields portable."""
    payload = json.dumps(data, separators=(",", ":"), allow_nan=False).encode() + b"\n"
    path.write_bytes(gzip.compress(payload, mtime=0))


def grain_palette(max_id):
    # Zero alone is dark; positive IDs never cycle back to the unassigned color.
    colors = ["#14271f"]
    for grain in range(1, max_id + 1):
        rgb = colorsys.hls_to_rgb((0.11 + grain * 0.61803398875) % 1,
                                 (0.53, 0.66, 0.59)[grain % 3], 0.48)
        colors.append("#" + "".join(f"{round(channel * 255):02x}" for channel in rgb))
    return colors


def read_vtk(path, slice_y):
    """Stream the full VTK; retain three exact planes and visible exterior faces."""
    result, dimensions, spacing = {}, None, None
    with path.open() as stream:
        for line in stream:
            words = line.split()
            if not words:
                continue
            if words[0] == "DIMENSIONS":
                dimensions = [int(n) - 1 for n in words[1:4]]
                nx, ny, nz = dimensions
                planes = sorted({0, slice_y, ny - 1})
                plane_set = set(planes)
            elif words[0] == "SPACING":
                spacing = [float(n) for n in words[1:4]]
                if len(set(spacing)) != 1:
                    raise ValueError("Viewer requires cubic cells")
            elif words[0] == "SCALARS":
                name = words[1]
                next(stream)  # LOOKUP_TABLE default
                values, top, side, assigned = [], [], [], 0
                convert = int if name in ("grain_id", "cell_state") else float
                for k in range(nz):
                    for j in range(ny):
                        for i in range(nx):
                            value = convert(next(stream).strip())
                            if not math.isfinite(value):
                                raise ValueError(f"Non-finite {name} in {path}")
                            if name == "cell_state":
                                continue
                            if j in plane_set:
                                values.append(value)
                            if name == "grain_id":
                                assigned += value > 0
                                if k == nz - 1:
                                    top.append(value)
                                if i == nx - 1:
                                    side.append(value)
                if name != "cell_state":
                    result[name] = values
                if name == "grain_id":
                    result["surfaces"] = dict(top=top, side=side)
                    result["assigned_fraction"] = assigned / math.prod(dimensions)
    if dimensions is None or spacing is None or not {"grain_id", "temperature_K"} <= result.keys():
        raise ValueError(f"Incomplete VTK: {path}")
    return dimensions, spacing[0], planes, result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=ROOT / "output/fem_replay")
    parser.add_argument("--output", type=Path, default=ROOT / "docs/assets/demo")
    parser.add_argument("--peak-time-ms", type=float, default=0.4)
    args = parser.parse_args()
    summary = tomllib.loads((args.input / "summary.toml").read_text())
    middle = "peak" if (args.input / "peak.vtk").exists() else "midpoint"
    middle_time = summary.get("snapshot_time_ms", args.peak_time_ms) if middle == "peak" else summary["time_ms"] / 2
    frames, geometry = [], None
    stages = [("initial", 0.0), (middle, middle_time)]
    for name, time in (("first_cycle", 2.0), ("second_peak", 2.4)):
        if (args.input / f"{name}.vtk").exists():
            stages.append((name, time))
    stages.append(("final", summary["time_ms"]))
    for name, time in stages:
        dims, spacing, planes, fields = read_vtk(args.input / f"{name}.vtk", summary["slice_y_index"] - 1)
        if geometry is not None and geometry != (dims, spacing):
            raise ValueError("Snapshot geometry differs")
        geometry = dims, spacing
        frames.append(dict(name=name, time_ms=time, **fields))
    palette = grain_palette(summary["final_grains"])
    data = dict(dimensions=dims, slice_indices=planes, spacing_mm=spacing, frames=frames, summary=summary, palette=palette,
                provenance="Actual CalculiX temperature replay through the maintained thesis-derived Julia CA, including the mushy-neighbour growth restriction; illustrative, uncalibrated case.")
    dest = args.output / "grains.json.gz"
    dest.parent.mkdir(parents=True, exist_ok=True)
    write_gzip_json(dest, data)
    dest.with_suffix("").unlink(missing_ok=True)  # Replaces the former uncompressed export.
    # The scientific PNG panels provide a compact no-JavaScript fallback.
    # The full cell-by-cell SVG remains available in the local run directory.
    (dest.parent / "fem-replay.svg").unlink(missing_ok=True)
    # Stream one numerical plane at a time, including validation; no full history
    # allocation. Temperatures remain in the local CSV and the saved VTK stages.
    history_path = args.input / "slice_history.csv"
    stride = summary.get("slice_stride", 1)
    with history_path.open("rb") as source:
        source_hash = hashlib.file_digest(source, "sha256").hexdigest()
    evolution = dict(dimensions=[dims[0], dims[2]], spacing_mm=spacing, palette=palette,
        y_um=summary["slice_y_um"], peak_time_ms=middle_time, dt_ms=summary["dt_ms"],
        slice_stride=stride, steps=summary["steps"],
        heating_intervals_ms=summary.get("heating_intervals_ms", [[0.0, middle_time]]),
        source_sha256=source_hash, numerical_data="slice-data.jsonl.gz",
        provenance="Recorded CA grain IDs at a fixed x-z plane; no spatial or temporal interpolation. Full cell temperatures remain in the simulation CSV and VTK outputs.")
    timeline, frame = [], None
    with (dest.parent / "slice-data.jsonl.gz").open("wb") as compressed:
      with gzip.GzipFile(fileobj=compressed, filename="", mode="wb", mtime=0) as numerical:
        numerical.write((json.dumps(evolution, separators=(",", ":")) + "\n").encode())
        def finish_frame(frame):
            if len(frame["grain_id"]) != dims[0] * dims[2]:
                raise ValueError("Incomplete slice frame")
            if not math.isclose(frame["time_ms"], frame["step"] * summary["dt_ms"], abs_tol=1e-8):
                raise ValueError("Slice time disagrees with simulation timestep")
            numerical.write((json.dumps(frame, separators=(",", ":"), allow_nan=False) + "\n").encode())
            timeline.append(dict(step=frame["step"], time_ms=frame["time_ms"],
                grains=len(set(frame["grain_id"]) - {0}),
                assigned_fraction=sum(g > 0 for g in frame["grain_id"]) / len(frame["grain_id"])))
        with history_path.open(newline="") as stream:
            for row in csv.DictReader(stream):
                step = int(row["step"])
                if frame is None or step != frame["step"]:
                    if frame is not None:
                        finish_frame(frame)
                    frame = dict(step=step, time_ms=float(row["time_ms"]), grain_id=[])
                expected = len(frame["grain_id"])
                if (int(row["x_index"]), int(row["z_index"])) != (expected % dims[0] + 1, expected // dims[0] + 1):
                    raise ValueError("Slice cells are not ordered x-fastest")
                if int(row["y_index"]) != summary["slice_y_index"]:
                    raise ValueError("Slice plane changed during export")
                grain_id = int(row["grain_id"])
                if not 0 <= grain_id < len(palette) or not math.isfinite(float(row["temperature_K"])):
                    raise ValueError("Invalid grain ID or cell temperature")
                frame["grain_id"].append(grain_id)
        if frame is not None:
            finish_frame(frame)
    expected_steps = sorted({*range(0, summary["steps"] + 1, stride), summary["steps"]})
    if len(timeline) != summary["slice_frames"] or [f["step"] for f in timeline] != expected_steps:
        raise ValueError("Slice steps disagree with recorded output stride")
    evolution["frames"] = timeline
    (dest.parent / "evolution.json").write_text(json.dumps(evolution, separators=(",", ":"), allow_nan=False) + "\n")
    (dest.parent / "slice-data.json.gz").unlink(missing_ok=True)

    print(f"Exported {len(frames)} snapshots, {math.prod(dims)} cells each → {dest}")
    print(f"Exported {len(timeline)} actual 2D steps at y={evolution['y_um']} µm")


if __name__ == "__main__":
    main()
