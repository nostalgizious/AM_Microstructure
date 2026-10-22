#!/usr/bin/env python3
"""Scientific figures and frame images from retained FEM data and actual CA output."""
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import json
import re
import tomllib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap, to_rgba
import numpy as np
from PIL import Image
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/assets/demo"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False,
    "axes.spines.right": False, "svg.fonttype": "path", "svg.hashsalt": "ammicrostructure",
    "savefig.facecolor": "white", "axes.labelcolor": "#26372e", "text.color": "#26372e"})


def save(fig, name):
    fig.savefig(OUT / f"{name}.svg", metadata={"Date": None}, bbox_inches="tight")
    fig.savefig(OUT / f"{name}.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def main():
    global OUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=OUT)
    parser.add_argument("--run", type=Path, default=ROOT / "output/fem_replay")
    args = parser.parse_args()
    OUT = args.data
    OUT.mkdir(parents=True, exist_ok=True)
    data = json.loads(gzip.decompress((OUT / "grains.json.gz").read_bytes()))
    evolution = json.loads((OUT / "evolution.json").read_text())
    palette = data["palette"]
    nx, ny, nz = data["dimensions"]
    dx_um = data["spacing_mm"] * 1000
    ns = len(data["slice_indices"])
    y_index = data["slice_indices"].index(data["summary"]["slice_y_index"] - 1)
    max_id = len(palette) - 1
    cmap = ListedColormap(palette)
    norm = BoundaryNorm(np.arange(-.5, max_id + 1.5), cmap.N)
    extent = (0, nx * dx_um, 0, nz * dx_um)

    # Cell-resolved static panels: each pixel is a real simulated cell.
    fig, axes = plt.subplots(2, len(data["frames"]), figsize=(14, 5.5), layout="constrained", sharex=True, sharey=True)
    labels = {"initial": "Initial substrate", "peak": "Pass 1: heated", "first_cycle": "Pass 1: cooled", "second_peak": "Pass 2: heated", "final": "Pass 2: cooled"}
    for column, frame in enumerate(data["frames"]):
        label = labels[frame["name"]]
        ids = np.array(frame["grain_id"]).reshape((nx, ns, nz), order="F")[:, y_index, :].T
        temp = np.array(frame["temperature_K"]).reshape((nx, ns, nz), order="F")[:, y_index, :].T
        grains = axes[0, column].imshow(ids, origin="lower", extent=extent, cmap=cmap, norm=norm, interpolation="nearest")
        thermal = axes[1, column].imshow(temp, origin="lower", extent=extent, cmap="inferno", vmin=1500, vmax=1800, interpolation="nearest")
        axes[0, column].set_title(f"{label}\nt = {frame['time_ms']:.2f} ms", fontsize=10)
        axes[1, column].set_xlabel("x (µm)")
        for row in range(2):
            axes[row, column].set_xticks([0, 16, 32, 48, 64])
            axes[row, column].set_yticks([0, 10, 20, 30, 40])
    axes[0, 0].set_ylabel("Grain field\nz (µm)")
    axes[1, 0].set_ylabel("Temperature field\nz (µm)")
    fig.colorbar(grains, ax=axes[0, :], ticks=range(0, max_id + 1, max(1, max_id // 8)),
                 label="Grain ID (0: unassigned)", shrink=.9)
    fig.colorbar(thermal, ax=axes[1, :], label="Temperature (K)", shrink=.9)
    fig.suptitle(f"FEM-driven heating and regrowth · x–z slice at y = {evolution['y_um']:.2f} µm", fontsize=12)
    save(fig, "fem-grain-evolution")

    # Actual .dat nodal values: retain total times across both solver steps.
    dataset_name = data["summary"].get("thermal_dataset", "data/two_pass_heating_demo")
    dataset = ROOT / dataset_name
    mesh = tomllib.loads((dataset / "mesh.toml").read_text())
    initial = tomllib.loads((dataset / "initial.toml").read_text())
    shape = tuple(len(mesh[a]) for a in ("x_m", "y_m", "z_m"))
    node_ids = np.array(mesh["node_ids"]).reshape(shape, order="F")
    center = tuple(n // 2 for n in shape)
    top_x = [int(np.argmin(np.abs(np.array(mesh["x_m"]) * 1e6 - x))) for x in (24, 40)]
    selected = [node_ids[i, center[1], -1] for i in top_x] + [node_ids[center], node_ids[center[0], center[1], 0]]
    times, blocks = [0.0], [dict(zip(initial["node_ids"], initial["temperature_K"]))]
    for line in (dataset / "thermal.dat").read_text().splitlines():
        match = re.search(r"temperatures for set NALL and time\s+(\S+)", line)
        if match:
            times.append(float(match[1].replace("D", "E")) * 1000)
            blocks.append({})
        elif len(times) > 1 and re.match(r"^\s*\d+\s+[-+0-9.]", line):
            node, value = line.split()
            blocks[-1][int(node)] = float(value.replace("D", "E"))
    assert len(times) == data["summary"]["solver_snapshots"] + 1 and all(len(b) == np.prod(shape) for b in blocks)
    with (args.run / "history.csv").open() as stream:
        history = list(csv.DictReader(stream))
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.8), layout="constrained")
    for node, label, color in zip(selected, ("Top at x = 24 µm", "Top at x = 40 µm", "Interior centre", "Bottom centre"), ("#a54028", "#915fa6", "#2c7381", "#47673e")):
        axes[0].plot(times, [b[int(node)] for b in blocks], label=label, color=color, lw=1.6)
    axes[0].axhline(1609, color="#667", ls="--", lw=.8, label="CA liquidus: 1609 K")
    axes[0].axhline(1533, color="#999", ls=":", lw=.8, label="CA solidus: 1533 K")
    axes[0].set(ylabel="Nodal temperature (K)", ylim=(1475, 1820), title="(a) Computed thermal history")
    axes[0].legend(fontsize=7, loc="upper center", bbox_to_anchor=(.5, -.20), frameon=False, ncol=2)
    axes[1].plot([float(r["time_ms"]) for r in history], [100 * float(r["assigned_fraction"]) for r in history], color="#316947", lw=1.8)
    axes[1].set(ylabel="Cells with an assigned grain ID (%)", ylim=(0, 105), title="(b) CA grain assignment")
    for ax in axes:
        for start, end in evolution["heating_intervals_ms"]:
            ax.axvspan(start, end, color="#f0d9bd", alpha=.5, zorder=-1)
            ax.axvline(end, color="#9b7852", lw=.8, ls="--")
        ax.set(xlim=(0, data["summary"]["time_ms"]), xlabel="Time (ms)")
        ax.grid(axis="y", alpha=.2)
    for cycle, (start, end) in enumerate(evolution["heating_intervals_ms"], 1):
        axes[1].text((start + end)/2, 8, f"Pass {cycle}", ha="center", fontsize=8)
    save(fig, "fem-thermal-history")

    # A static three-dimensional figure remains usable without JavaScript.
    ids = np.array(data["frames"][-1]["grain_id"]).reshape((nx, ns, nz), order="F")
    fig = plt.figure(figsize=(7, 5.2))
    ax = fig.add_subplot(projection="3d")
    # Only the three visible exterior faces: no internal voxel geometry or cell
    # gridlines that would obscure grain boundaries on the refined CA grid.
    top = np.array(data["frames"][-1]["surfaces"]["top"]).reshape((nx, ny), order="F")
    side = np.array(data["frames"][-1]["surfaces"]["side"]).reshape((ny, nz), order="F")
    quads, colors = [], []
    for x in range(nx):
        for y in range(ny):
            quads.append([(x,y,nz), (x+1,y,nz), (x+1,y+1,nz), (x,y+1,nz)])
            colors.append(to_rgba(palette[top[x,y]]))
    for x in range(nx):
        for z in range(nz):
            quads.append([(x,0,z), (x+1,0,z), (x+1,0,z+1), (x,0,z+1)])
            colors.append(to_rgba(palette[ids[x,0,z]]))
    for y in range(ny):
        for z in range(nz):
            quads.append([(nx,y,z), (nx,y+1,z), (nx,y+1,z+1), (nx,y,z+1)])
            color = np.array(to_rgba(palette[side[y,z]])); color[:3] *= .8
            colors.append(color)
    ax.add_collection3d(Poly3DCollection(np.array(quads) * dx_um, facecolors=colors,
                                        linewidths=0, edgecolors="none", antialiased=False, rasterized=True))
    ax.set(xlim=(0,nx*dx_um), ylim=(0,ny*dx_um), zlim=(0,nz*dx_um))
    ax.set(title=f"Final grain field after two heating passes · t = {data['summary']['time_ms']:.0f} ms")
    ax.set_box_aspect((nx, ny, nz)); ax.view_init(elev=23, azim=-60)
    ax.set_xticks([0, 32, 64]); ax.set_yticks([0, 12]); ax.set_zticks([0, 20, 40])
    fig.subplots_adjust(left=.03, right=.85, bottom=.12, top=.90)
    # Explicit figure labels avoid clipped 3D-axis labels in SVG renderers.
    ax.text2D(.48, .02, "x (µm)", ha="center", transform=ax.transAxes)
    ax.text2D(.94, -.02, "y (µm)", ha="center", transform=ax.transAxes)
    ax.text2D(1.0, .55, "z (µm)", ha="center", rotation=90, transform=ax.transAxes)
    save(fig, "grain-volume")

    # One downloadable image for every real CA output step, with no blending.
    frame_dir = OUT / "evolution"
    frame_dir.mkdir(exist_ok=True)
    rgba = np.array([[int(color[i:i+2], 16) for i in (1, 3, 5)] + [255]
                     for color in palette], dtype=np.uint8)
    for old_frame in frame_dir.glob("frame-*.png"):
        old_frame.unlink()
    magnification = max(1, 1024 // nx)
    with gzip.open(OUT / "slice-data.jsonl.gz", "rt") as stream:
        header = json.loads(next(stream))
        assert header["source_sha256"] == evolution["source_sha256"]
        for line in stream:
            frame = json.loads(line)
            ids = np.array(frame["grain_id"]).reshape((nz, nx))
            pixels = rgba[np.flipud(ids)]
            im = Image.fromarray(pixels).resize((nx * magnification, nz * magnification), Image.Resampling.NEAREST)
            im.save(frame_dir / f"frame-{frame['step']:03d}.png", optimize=True)
            if abs(frame["time_ms"] - 2.0) < 1e-8:
                im.save(OUT / "first-cycle.png", optimize=True)
    files = sorted([*OUT.glob("*.svg"), *OUT.glob("*.png"), *frame_dir.glob("*.png"), *OUT.glob("*.gz"), OUT / "evolution.json"])
    manifest = dict(generator="tools/make_demo_figures.py", dataset=dataset_name,
        source_csv_sha256=evolution["source_sha256"], slice_y_um=evolution["y_um"],
        ca_dimensions=data["dimensions"], spacing_um=dx_um,
        initial_grains=data["summary"]["initial_grains"], final_grains=data["summary"]["final_grains"],
        growth_rule=data["summary"].get("growth_rule", "not recorded"),
        rng_contract=data["summary"].get("rng_contract", "not recorded"),
        worker_threads=data["summary"].get("worker_threads"),
        dt_ms=evolution["dt_ms"],
        frame_count=len(evolution["frames"]), units=dict(length="µm", time="ms", temperature="K"),
        scope="Actual two-pass demonstration outputs; categorical grain IDs, uncalibrated thermal specimen; no temporal image interpolation.",
        files=[dict(file=p.relative_to(OUT).as_posix(), sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files])
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Generated 3 figure pairs (SVG/PNG) and {len(evolution['frames'])} actual slice images.")


if __name__ == "__main__":
    main()
