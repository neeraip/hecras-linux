# HEC-RAS 2D on Linux

Run HEC-RAS 2D unsteady simulations on Linux or in Docker, without Windows or
the HEC-RAS GUI, using HEC's official Linux engines (`RasGeomPreprocess`,
`RasUnsteady`). Engines for **HEC-RAS 6.6**, **7.0.1** and **7.0.2** are included.

`ras_preprocess.py` prepares a project for the engines. What it has to do
depends on what the project already contains, which gives three workflows:

| | Input | What the preprocessor does | Status |
|---|---|---|---|
| **A** | GUI-computed project (`g01.hdf` with hydraulic tables, and `p01.hdf`) | Copies the geometry and strips the old results | **Validated.** Matches the Windows GUI to within 0.0002 ft on every example, with all three engines |
| **B** | RASMapper mesh in `g01.hdf`, geometry not yet computed | Computes the property tables as RASMapper does (`property_tables.py`) | **Tables match RASMapper's** value for value on the examples (see [accuracy](#workflow-b-and-c-accuracy)) |
| **C** | Only the `.g01` text file (seed points, perimeter, BC lines) | Builds the mesh, then does B | **Experimental.** Mesh close to RASMapper's but not identical, so results differ |

The workflow is picked automatically from the files present.

## Quick start

```bash
git clone https://github.com/neeraip/hecras-v66-linux.git
cd hecras-v66-linux                 # engines arrive via git LFS (~900 MB)
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
RASMapper computed. They do not test the Python code in Workflows B and C.

## Workflow B and C accuracy

`property_tables.py` is a port of RASMapper's "Compute Property Tables" code
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
  row of the infiltration layer's variables.

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
libraries; none appears in the examples): *Multiple Face Mann n* (refused with
an error), cell-averaged land cover (*Composite LC*), calibration and override
regions, several terrain rasters with priorities, terrain modifications, and
2D structures and connections. One 2D flow area and US customary units only.

Workflow C builds its own mesh, which is close to RASMapper's but not
identical: from the same seed points and perimeter it produces the same
interior cells and within 0–4 faces of RASMapper's count (the difference is in
how a few perimeter corners are split). Its tables are computed exactly as in
B, but on that mesh, so its results differ from the GUI's where the meshes do.

### Limitations of Workflows B and C

- **Boundary data.** The engines take hydrographs from the DSS file named in
  the plan. Without it, pass `--event-conditions <HDF>` to copy boundary and
  initial conditions from an existing plan HDF (for example the GUI's
  `.p01.hdf`); their face lists are rebuilt for the new mesh, with whole faces
  only (RASMapper can use a fraction of the last face of a BC line).
- **Mesh (C).** Uniform seed points from the `.g01` only: no breaklines,
  refinement regions or structures. Cells with more than 8 faces are reported
  but not repaired.
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
`p01.tmp.hdf` from the `.p01` text. *Multiple Face Mann n* is not supported.

### Workflow C — no g01.hdf

**When:** there is no `g01.hdf`. The preprocessor reads the seed points,
perimeter and BC lines from the `.g01` text, clips each seed's Voronoi cell to
the perimeter, and stitches the cells into a mesh that follows RASMapper's
conventions (perimeter faces between face points, one ghost cell per perimeter
face, counter-clockwise cells). It then continues as Workflow B.

## Repository layout

```
ras_preprocess.py          preprocessor (Workflows A / B / C)
property_tables.py         RASMapper's 2D property-table computation, ported
hecras_run.py              runs the engines natively or in Docker, checks the result
run_simulation.sh          zip in → preprocess → simulate → zip out
tools/
  make_workflow_input.py   derive a Workflow B or C input from a GUI project
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
every example (`make test-fast` checks this, about 4 minutes). Workflow B and C
tests check that the plan runs; their accuracy tests are marked as expected
failures until they are measured on every example. The run ends with a
table of volume error and WSE difference for every run.

## Adding another engine version

Put the three executables in `bin/<version>/` and compare its libraries with
`libs/`. So far every release has kept the existing libraries unchanged and at
most added new ones (7.0.2 added four MKL VML libraries), which can simply be
copied into `libs/`. If a release ever changes an existing library, give it
its own `libs/` and point `hecras_run.py` at it. Then
`pytest --hecras-version <version>` and
`docker build --build-arg HECRAS_VERSION=<version> -t hecras-linux:<version> .`.
