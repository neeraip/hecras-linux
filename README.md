# hecras-linux

Run HEC-RAS 2D unsteady simulations on Linux or in Docker, without Windows or
the HEC-RAS GUI, using HEC's official Linux engines (`RasGeomPreprocess`,
`RasUnsteady`). Engines for **HEC-RAS 6.6**, **7.0.1** and **7.0.2** are included.

The one step those engines can't do is RASMapper's "Compute Property Tables".
`hecras_linux.property_tables` is a port of it that reproduces RASMapper's
tables byte for byte, so a mesh built in RASMapper, or by any tool that writes
one, can run on Linux.

`ras_preprocess.py` (`hecras-preprocess` when installed) prepares a project
for the engines. What it has to do depends on what the project already
contains, which gives two workflows:

| | Input | What the preprocessor does | Status |
|---|---|---|---|
| **A** | GUI-computed project (`g01.hdf` with hydraulic tables, and `p01.hdf`) | Copies the geometry and strips the old results | **Validated.** Matches the Windows GUI to within 0.0002 ft on every example, with all three engines |
| **B** | A mesh in `g01.hdf`, geometry not yet computed | Computes the property tables as RASMapper does (`property_tables.py`) | **Tables match RASMapper's** byte for byte on the examples (see [accuracy](#workflow-b-accuracy)) |

The workflow is picked automatically from the files present. Building the mesh
itself (from a perimeter, cell size and breaklines) is out of scope: use
RASMapper, or a mesh generator that writes the geometry HDF.

## Quick start

```bash
git clone https://github.com/neeraip/hecras-linux.git
cd hecras-linux                     # engines arrive via git LFS (~900 MB)
pip install -r requirements.txt     # GDAL: see requirements.txt
make test-a                         # Workflow A on every example
```

Run one project from a zip:

```bash
./run_simulation.sh MyProject.zip                       # → MyProject_results.zip
HECRAS_VERSION=7.0.2 ./run_simulation.sh MyProject.zip  # use the 7.0.2 engines
```

Or step by step:

```bash
python3 ras_preprocess.py Examples/BEC --output-dir /tmp/bec
cp Examples/BEC/BEC.b01 Examples/BEC/BEC.x01 /tmp/bec/
python3 hecras_run.py /tmp/bec BEC --version 6.6        # writes /tmp/bec/BEC.p01.hdf
python3 tools/compare_results.py /tmp/bec/BEC.p01.hdf Examples/BEC/BEC.p01.hdf
```

### Using the property tables from another tool

```bash
pip install "hecras-linux @ git+https://github.com/neeraip/hecras-linux@<tag>"
```

```python
from hecras_linux import property_tables, preprocess

terrain = property_tables.TerrainLayer("Terrain/Terrain.hdf")    # or a .tif
s, z = property_tables.face_profile(terrain, face_polyline, tol)
rows = property_tables.face_table(s, z, mannings_n, area_tol, conveyance_ratio, laminar_depth)
table = property_tables.cell_table(terrain, cell_polygon, area, lowest_face, vol_tol, min_area_fraction)
```

`preprocess` does the whole of Workflow B on a project folder, as the
command line does. Installing the package doesn't bring the engines; set
`HECRAS_HOME` to a clone (with `bin/` and `libs/`) to run them through
`hecras_linux.run`.

### Running the engines: native or Docker

`hecras_run.py` runs the engines **natively** on Linux x86_64 (tested on
Ubuntu 24.04; glibc 2.28 or newer). On macOS, Windows or ARM use **Docker**:

```bash
make docker                                   # builds hecras-linux:6.6, :7.0.1 and :7.0.2
HECRAS_RUNNER=docker ./run_simulation.sh MyProject.zip
python3 hecras_run.py /tmp/bec BEC --runner docker
```

The engines exit with status 0 even when they stop on an input error, so
`hecras_run.py` checks the log and only reports success when the simulation
finished and wrote results.

### Choosing an engine version

The 7.0.1 and 7.0.2 engines read 6.6 projects unchanged and give the same
results on every example (7.0.1 and 7.0.2 are identical to each other). Plans
they write can't be opened by HEC-RAS 6.6, so the default stays 6.6. Use 7.0.2
when the results will be opened in HEC-RAS 7.

## Validated results (Workflow A)

Linux results against the Windows GUI's own results for the same plan:

| Example | Cells | Volume error | Max WSE diff, 6.6 | 7.0.1 | 7.0.2 |
|---|---|---|---|---|---|
| AFTER_RUN | 375 | 0.099 % | 0.0000 ft | 0.0000 ft | 0.0000 ft |
| BEC | 2,315 | 0.002 % | 0.0001 ft | 0.0001 ft | 0.0001 ft |
| BEC_WO_Infiltration | 2,622 | 0.001 % | 0.0000 ft | 0.0000 ft | 0.0000 ft |
| VA | 1,221 | 0.00001 % | 0.0000 ft | 0.0000 ft | 0.0000 ft |
| Muncie | 3,146 | 0.049 % | 0.0001 ft | 0.0002 ft | 0.0002 ft |

These show that the Linux engines reproduce Windows HEC-RAS on geometry that
RASMapper computed. They do not test the Python code in Workflow B.

## Workflow B accuracy

`hecras_linux/property_tables.py` is a port of RASMapper's "Compute Property Tables" code
(RasMapperLib, HEC-RAS 7.0.2; HEC-RAS is public-domain USACE software). It
follows RASMapper step by step, down to where RASMapper rounds to 32-bit
floats:

- **Terrain** is read the way RASMapper reads its terrain GeoTIFF: in "render
  blocks" between four pixel centres, each split into four triangles meeting
  at the block centre (the mean of the four corners), located through
  RASMapper's tile grid. This is not bilinear interpolation, and the
  difference shows in every table. Blocks flagged in the terrain HDF's masks
  as bordering no data are bridged with RASMapper's stitch TIN.
- **Faces:** the profile is sampled where the face crosses render-block edges
  and diagonals, thinned per terrain tile and again overall (Douglas–Peucker,
  *Face Profile Tolerance*). The area–elevation table is computed at the
  profile's elevations, refined and thinned on conveyance (*Face Area
  Conveyance Ratio*, *Laminar Depth*). The low-elevation centroid is the
  centroid of the bottom 5 % of the face's area.
- **Cells:** each cell is resampled onto its own raster of up to 10,000 pixels,
  binned into a 256-bin area histogram, integrated to volume and thinned
  (*Cell Volume Tolerance*). The bottom is raised to where *Cell Minimum Area
  Fraction* of the cell is wet, and lowered to the cell's lowest face.
- **Land cover and infiltration:** Manning's *n* and classes at each cell
  centre and face midpoint, off-raster points taking the NoData class.
  Infiltration classes combine the land cover and soils classes, and select a
  row of the infiltration layer's variables. With *Composite LC*, a cell's
  *n* is the average over its raster instead. With *Multiple Face Mann n*,
  each face's *n* follows the land cover along it: the face is split where
  *n* changes, and each table row's *n* combines the wetted parts. Manning's
  *n* calibration regions override *n* by land cover class inside their
  polygons, at points and along faces.

Workflow B is tested by stripping the GUI's tables from each example
(`tools/make_workflow_input.py`) and rebuilding them.
`test_workflow_b_tables_match_rasmapper` compares every dataset RASMapper
writes for the 2D area with the GUI's: on all five examples they are
identical, byte for byte.

Engine results under Workflow B, 6.6 engine, GUI results as reference:

| Example | Volume error (GUI) | Max WSE difference |
|---|---|---|
| BEC | 0.0020 % (0.0020 %) | 0.000 ft on every cell |
| BEC_WO_Infiltration | 0.0014 % (0.0014 %) | 0.000 ft on every cell |

AFTER_RUN, VA and Muncie have not been re-run through the engine since the port.

Not ported yet (RASMapper's code for each is in RasMapperLib or its companion
libraries): calibration regions together with *Composite LC*, and regions
with several parts (both refused with an error); other override regions;
per-face variable averages under *Composite LC*; terrain modifications other
than polygons with one elevation, and faces crossing any modification (both
refused with an error); 2D structures and connections. One 2D flow area and US
customary units only. Terrains with several rasters are read in RASMapper's
priority order, with its stitch TIN across their seams, and a terrain
modification polygon with one elevation (set, add, take higher or take lower)
is applied to the cells under it.

Beyond the shipped examples, every table of HEC's Weise_2D example and of
BaldEagleCrkMulti2D g09 (two rasters; tables recomputed in RASMapper 7.0.2)
comes out bit-exact, and so does every table of Chippewa_2D (*Composite LC*
and *Multiple Face Mann n*, tables recomputed in RASMapper 7.0.2). On a
production model (Multiple Face Mann n, a calibration region and a terrain
modification), every table compared is identical: all 112,246 face profiles,
area-elevation tables, minimum elevations and low-elevation centroids, all
47,103 cell volume tables and Manning's *n* at all 48,006 cell centres.

[`docs/rasmapper-coverage.md`](docs/rasmapper-coverage.md) measures coverage
against RASMapper itself: every method its table computation can reach, read
from the IL of HEC's DLLs, with its status and our counterpart.

### Limitations of Workflow B

- **Boundary data.** The engines take hydrographs from the DSS file named in
  the plan. Without it, pass `--event-conditions <HDF>` to copy boundary and
  initial conditions from an existing plan HDF (for example the GUI's
  `.p01.hdf`); their face lists are rebuilt for the new mesh, with whole faces
  only (RASMapper can use a fraction of the last face of a BC line).
- **One 2D flow area**, US Customary units, unsteady flow only.

## How the workflows work

### Workflow A — GUI-computed project

**When:** `g01.hdf` already has hydraulic tables **and** a `p01.hdf` exists.
The preprocessor copies `g01.hdf` and writes `p01.tmp.hdf` from `p01.hdf`
without its `Results` group. `RasGeomPreprocess` and `RasUnsteady` then
recompute from scratch.

### Workflow B — RASMapper mesh, tables not computed

**When:** `g01.hdf` exists but has no hydraulic tables. The preprocessor reads
the mesh, computes the property tables with RASMapper's algorithm
(`property_tables.py`, settings from the 2D area's attributes), assigns
Manning's *n*, percent impervious and infiltration from the land cover and
soils layers, assigns perimeter faces to BC lines, and assembles
`p01.tmp.hdf` from the `.p01` text.

## Repository layout

```
hecras_linux/
  property_tables.py       RASMapper's 2D property-table computation, ported
  preprocess.py            preprocessor (Workflows A and B)
  run.py                   runs the engines natively or in Docker, checks the result
ras_preprocess.py          command line for hecras_linux.preprocess
hecras_run.py              command line for hecras_linux.run
pyproject.toml             the hecras-linux package
run_simulation.sh          zip in → preprocess → simulate → zip out
tools/
  make_workflow_input.py   derive a Workflow B input from a GUI project
  rasmapper_callgraph/     every RASMapper method the table computation reaches (IL)
  rasmapper_coverage.py    classify them: ported, not ported, ... → docs/rasmapper-coverage.md
  compare_results.py       compare max WSE between two plan HDFs
tests/test_pipeline.py     pytest: all workflows, all examples, all engines
Examples/                  six HEC-RAS 6.6 projects (see Examples/README.md)
bin/6.6/ 7.0.1/ 7.0.2/     HEC-RAS Linux engines (git LFS)
libs/                      shared libraries for all engine versions (git LFS)
Dockerfile                 hecras-linux:<version> image
RAS_v.6.6_Linux.pdf        HEC's documentation for the Linux engines
remove_HDF5_Results_Sed.py strip results from a plan HDF (as Workflow A does)
```

## Tests

```bash
make test-fast                         # preprocessing only, ~5 min
make test                              # all workflows, every engine under bin/
pytest --hecras-version 7.0.2          # one engine
pytest --runner docker                 # through the Docker images
```

Workflow A tests must pass. Workflow B's tables must equal RASMapper's on
every example (`make test-fast` checks this, about 4 minutes). Workflow B's
engine tests check that the plan runs; its accuracy tests are marked as
expected failures until they are measured on every example. The run ends with a
table of volume error and WSE difference for every run.

## Adding another engine version

Put the three executables in `bin/<version>/` and compare its libraries with
`libs/`. So far every release has kept the existing libraries unchanged and at
most added new ones (7.0.2 added four MKL VML libraries), which can simply be
copied into `libs/`. If a release ever changes an existing library, give it
its own `libs/` and point `hecras_run.py` at it. Then
`pytest --hecras-version <version>` and
`docker build --build-arg HECRAS_VERSION=<version> -t hecras-linux:<version> .`.
