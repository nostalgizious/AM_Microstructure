# Demo engineering notes

This file collects implementation, runtime, provenance, and maintenance details for the demonstrations. The public guides stay focused on the cellular automaton and its thermal inputs.

## Numerical conventions and implementation

### Nucleation

The shared kernel integrates a Gaussian over the change in thermal undercooling, `T_liquidus − T`. The integral is multiplied by cell volume `dx³`; the historical normalization omits the conventional `1 / sigma` factor. The configuration's `N0` therefore acts as a model parameter tied to this normalization. Nucleation probability is capped at 0.005 per substep. Gaussian integration uses bounded composite Simpson quadrature, truncating tails at 12 standard deviations.

These conventions document the current implementation. They do not establish calibrated nucleation parameters. Research context and citations are in the [research guide](../docs/research.md).

### Thermal demonstrations

The directional-solidification example prescribes

```text
T(z, t) = T_liquidus + superheat + gradient * z - cooling_rate * t
z = (k - 1/2) * dx
```

The bottom layer remains 5 K below solidus and starts with seeded grains. The 32 × 6 × 24 grid spans 0.064 × 0.012 × 0.048 mm. Distances use mm, temperatures K, elapsed time ms, and orientation angles degrees. The growth-velocity expression has an SI velocity scale; a numerical velocity in m/s has the same value in mm/ms. The thermal field is prescribed rather than solved by a heat equation.

The default two-pass FEM dataset represents a single existing solid specimen, 0.064 × 0.012 × 0.040 mm, with 256 DC3D8 elements and 459 nodes. Two stationary Gaussian top-temperature profiles have σ = 10 µm and centres at x = 24 and 40 µm. Each rises from 1500 K to 1800 K over 0.2 ms, holds to 0.4 ms, then cools to 2 ms; the second cycle runs from 2 to 4 ms. The bottom is held at 1500 K. During cooling, a top film uses a 100,000 W/(m² K) coefficient and a 1500 K environment. The material parameters are conductivity 25 W/(m K), density 7800 kg/m³, and heat capacity 600 J/(kg K).

The retained solver history has 200 snapshots at 0.02 ms. The showcase CA preset uses a 384 × 72 × 240 grid, 1/6 µm spacing, 64 initial substrate grains, and 1,600 steps of 0.0025 ms. The `quick` and `fine` presets use coarser CA grids over the same saved thermal solution. Prescribed temperature boundaries do not define laser energy or power. The case does not include element activation, layer deposition, phase-change heat, mechanics, radiation, or melt flow.

The one-way handoff stores mesh axes in metres, node temperatures in kelvin, and time in seconds. The adapter maps cell centres by trilinear spatial interpolation and interpolates between saved time points; it converts lengths to mm and time to ms before calling the CA. Replay rejects missing mesh nodes, cells outside the structured mesh, and times outside the saved interval.

## Kernel restoration and state management

The maintained implementation retains the shared growth kernel and wraps it with configuration, field checks, persistent state, and output handling. The restoration work addressed the following geometry and state issues:

- Cell centres are initialized from voxel geometry before growth.
- Capture checks the complete rotated octahedron, not a single face half-space.
- Inverse Z–X–Z rotation reverses the order of the angles as well as their signs.
- Closest-corner selection uses full spatial distances.
- Candidate traversal covers the second Moore neighbourhood and continues after a failed candidate.
- Captures are committed after candidate traversal, preventing within-sweep chain reactions.
- Rapid cooling below solidus can clear stranded liquid labels.
- Octahedral growth requires a mushy-state neighbour in the first Moore neighbourhood (the 26 adjacent cells, clipped at domain boundaries).
- Configuration construction, field validation, and refresh preserve the intended parameters, including `multilayer`.

The second-Moore capture search and first-Moore growth eligibility rule serve separate purposes. The persistent `SimulationState` retains IDs, states, orientations, octahedra, prior temperature, time, and its random generator. Nucleation draws and grain IDs follow canonical cell order; capture uses deterministic per-target random streams and deferred state commits. One- and eight-worker results are compared within the same Julia version. Historical task-RNG outputs can differ.

`step!(state, Tnext; dt_ms, skip_cold=true)` can skip microstructure updates when all cells are assigned solid and both interval endpoints remain strictly below solidus. The prior temperature and simulation time still advance. Set `skip_cold=false` to run the complete kernel. Adapters should split at intermediate thermal extrema because endpoint temperatures cannot reveal an unsampled crossing.

## Saved-data provenance

Each dataset directory contains its solver deck, raw thermal output, explicit mesh, initial field, CA configuration, solver log, and provenance metadata. The provenance files record source identities and SHA256 checksums; replay checks required dataset files before advancing. Raw thermal output is consumed with the explicit initial snapshot. The CA state is reconstructed from its saved configuration and seed rather than serialized Julia objects.

The `reference_run.toml` files record example summaries. They are specific to their configuration and runtime. Dataset directories and their readmes are linked from the [FEM guide](../docs/fem-integration.md).

## CalculiX runtime setup

The optional fresh-solve path needs CalculiX CrunchiX (`ccx`); the saved-data replay path needs only Julia. On Debian or Ubuntu, the solver package is commonly installed with:

```sh
sudo apt-get update
sudo apt-get install calculix-ccx
```

Select the executable explicitly when needed:

```sh
julia --project=. --threads=8,0 integrations/calculix/run.jl \
  --ccx /absolute/path/to/ccx --output output/calculix
```

A Fedora run using a user-supplied ccx 2.23 binary required `libgfortran.so.4`. A compatible runtime was unpacked from a user-provided Ubuntu `libgfortran4_7.5.0-3ubuntu1~18.04_amd64.deb` into `/tmp/am-ccx-runtime`; the solve then ran with:

```sh
LD_LIBRARY_PATH=/tmp/am-ccx-runtime/usr/lib/x86_64-linux-gnu \
CCX=/home/ferdinand/Downloads/CalculiX/ccx_2.23/src/ccx_2.23 \
  julia --project=. --threads=8,0 integrations/calculix/run.jl
```

This was a machine-specific runtime workaround. The solver binary and extracted runtime are not included in the repository. A failed solver invocation is reported as an error; no generated thermal history is substituted.

### Mesh loading and fracture analysis

The tetrahedral mesh reader, barycentric mapping and grain-initialization routines are in [`src/coupling/`](../src/coupling/). The HDF5 mesh loaders require HDF5.jl and use the original host array layouts. The host solver and its C callback registration are not included. The supported thermal handoff is described in the [FEM guide](../docs/fem-integration.md).

The thesis fracture-parameter and plotting routines are in [`analysis/extractFractureParameters.jl`](../analysis/extractFractureParameters.jl). Their optional dependencies are HDF5, StatsBase, Plots, LaTeXStrings and LsqFit. The HDF5 routines retain the original study filenames and directory settings; set these for the data being analysed before calling them.

## Figure and website maintenance

The research homepage is static HTML; MkDocs builds the documentation using [`mkdocs.yml`](../mkdocs.yml). Python is needed to rebuild it; Julia and CalculiX are not required for a documentation-only build. The following commands were used for local maintenance:

```sh
python3 -m venv .venv-docs
.venv-docs/bin/pip install -r requirements-docs.txt
.venv-docs/bin/python tools/build_site.py
python3 tools/check_site.py
python3 -m http.server 8000 --directory website/dist
```

The complete build combines the homepage with MkDocs navigation, search and documentation. To preview only the documentation, use `mkdocs serve`; the full-site preview above also serves the homepage and source download.

Figure generation uses simulation output and optional plotting dependencies:

```sh
julia --project=. --threads=8,0 examples/fem_replay/run.jl
python3 tools/export_demo_data.py
python3 -m venv .venv-figures
.venv-figures/bin/pip install -r requirements-figures.txt
.venv-figures/bin/python tools/make_demo_figures.py
.venv-figures/bin/python tools/check_demo_frames.py
```

This produces numerical SVG/PNG figures, saved playback frames, and a manifest under `docs/assets/demo/`. The frame checker compares displayed pixels with exported grain IDs. Figures and playback manifests are linked from [Figures](../docs/figures.md).

Source-document figure extraction uses Poppler's `pdftoppm` and local PDFs in the ignored `references/local/` directory:

```sh
python3 tools/extract_figures.py
```

The extraction manifest at [`docs/assets/research/manifest.json`](../docs/assets/research/manifest.json) records pages, crop coordinates, attribution, rights, and hashes. The source filenames are listed in [`references/README.md`](../references/README.md).

The checked-in GitHub Actions workflows cover Julia regression tests, optional stock-CalculiX integration, static-site checks, and manual Pages deployment. The deployment workflow publishes the combined `website/dist/` artifact on pushes to `main` and on manual runs. Documentation builds use MkDocs strict mode and the complete site is checked for broken local links.

### GitHub Pages

From a committed `AM_Microstructure` checkout, run `bash tools/publish_github.sh` with an authenticated GitHub CLI. This creates the public repository, pushes `main`, enables GitHub Pages, sets the repository website URL and starts deployment. Existing `origin` repositories must also be named `AM_Microstructure`.

For manual setup:

1. Push the repository, including `.github/workflows/`, the saved datasets and figure assets.
2. In the repository's **Settings → Pages**, select **GitHub Actions** as the source.
3. Push to `main` or run **Actions → Deploy showcase → Run workflow**. The deployment job reports the published URL.
4. Add that URL to the repository's **About → Website** field.

The build uses relative internal URLs for deployment at a domain root or a project path. For another static host, upload the contents of `website/dist/`. The downloadable source archive explicitly includes the thesis PDF and excludes generated simulation output, local journal PDFs, runtime packages, environments and Git metadata.

## Recorded verification snapshot

The following machine results are retained as engineering history, not as general performance guarantees. On Fedora with Julia 1.13.2 and eight worker threads, the corrected two-pass showcase advanced 6,635,520 CA cells through 1,600 steps. One recorded run took 297.94 s and reached 2,237,372 KiB maximum resident memory. It skipped 1,188 cold steps, performed 1,890,432 octahedral captures, and ended with 64 grains and 100% assigned cells. It produced 401 recorded x–z sections over 4 ms. Assignment is a grain-label statistic.

The regression suite at that point reported 299 assertions, including geometry, melting, quadrature, seeded initialization, field checks, interpolation, units, schedules, and cold intervals followed by remelting. A separate saved-data replay produced matching array hashes and next random draws with one and eight workers on that Julia version. These records describe the software and data checks made at the time; they are not physical validation claims.
