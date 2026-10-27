# HEC-RAS 2D on Linux

Run HEC-RAS 2D unsteady simulations on Linux or in Docker, without Windows or
the HEC-RAS GUI, using HEC's official Linux engines (`RasGeomPreprocess`,
`RasUnsteady`). Engines for **HEC-RAS 6.6** and **HEC-RAS 7.0.1** are included.

`ras_preprocess.py` prepares a project for the engines. What it has to do
depends on what the project already contains, which gives three workflows:

| | Input | What the preprocessor does | Status |
|---|---|---|---|
| **A** | GUI-computed project (`g01.hdf` with hydraulic tables, and `p01.hdf`) | Copies the geometry and strips the old results | **Validated.** Matches the Windows GUI to within 0.0002 ft on every example, with both engines |
| **B** | RASMapper mesh in `g01.hdf`, geometry not yet computed | Computes the hydraulic property tables in Python | **Experimental.** Runs, but results differ from the GUI (see [accuracy](#workflow-b-and-c-accuracy)) |
| **C** | Only the `.g01` text file (seed points, perimeter, BC lines) | Builds the mesh, then does B | **Experimental.** Mesh matches RASMapper's; tables have B's limitations |

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
HECRAS_VERSION=7.0.1 ./run_simulation.sh MyProject.zip  # use the 7.0.1 engines
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
make docker                                   # builds hecras-linux:6.6 and hecras-linux:7.0.1
HECRAS_RUNNER=docker ./run_simulation.sh MyProject.zip
python3 hecras_run.py /tmp/bec BEC --runner docker
```

The engines exit with status 0 even when they stop on an input error, so
`hecras_run.py` checks the log and only reports success when the simulation
finished and wrote results.

### Choosing an engine version

The 7.0.1 engines read 6.6 projects unchanged and give the same results on
every example. Plans they write can't be opened by HEC-RAS 6.6, so the default
stays 6.6. Use 7.0.1 when the results will be opened in HEC-RAS 7.

## Validated results (Workflow A)

Linux results against the Windows GUI's own results for the same plan:

| Example | Cells | Volume error | Max WSE diff, 6.6 | Max WSE diff, 7.0.1 |
|---|---|---|---|---|
| AFTER_RUN | 375 | 0.099 % | 0.0000 ft | 0.0000 ft |
| BEC | 2,315 | 0.002 % | 0.0001 ft | 0.0001 ft |
| BEC_WO_Infiltration | 2,622 | 0.001 % | 0.0000 ft | 0.0000 ft |
| VA | 1,221 | 0.00001 % | 0.0000 ft | 0.0000 ft |
| Muncie | 3,146 | 0.049 % | 0.0001 ft | 0.0002 ft |

These show that the Linux engines reproduce Windows HEC-RAS on geometry that
RASMapper computed. They do not test the Python code in Workflows B and C.

## Workflow B and C accuracy

Workflow B is tested by stripping the GUI's tables from each example
(`tools/make_workflow_input.py`) and letting Python rebuild them; Workflow C by
deleting the `g01.hdf` altogether. Both runs use the GUI's boundary data, so
any difference comes from the geometry.

Measured with the 6.6 engine (7.0.1 agrees to within 0.03 ft); WSE
differences are over all interior cells, GUI results as reference:

| Example | Workflow B: volume error | B: mean / worst WSE diff | Workflow C: volume error | C: mean / worst WSE diff |
|---|---|---|---|---|
| AFTER_RUN | 0.06 % | 0.12 / 3.3 ft | 0.06 % | 0.12 / 3.3 ft |
| BEC | 0.30 % | 0.91 / 17.8 ft | 0.30 % | 0.91 / 17.8 ft |
| BEC_WO_Infiltration | 0.34 % | 1.33 / 23.9 ft | 0.36 % | 1.32 / 23.5 ft |
| VA | 0.13 % | 2.41 / 31.2 ft | **29.7 %** | 2.57 / 31.2 ft |
| Muncie | **197.6 %** | 0.13 / 7.4 ft | **197.4 %** | 0.13 / 7.4 ft |

For comparison, Workflow A on the same examples: under 0.1 % volume error and
at most 0.0002 ft. B and C agree closely with each other, so the error comes
from the shared property-table code rather than from the C mesh. Two cases
stand out. Muncie's volume accounting breaks down in both; it is the only
example with both infiltration and precipitation, which is where to look
first. VA under C loses the partial face fractions RASMapper uses for its
boundary line, which wraps round the whole area (B keeps the GUI's).

Why B differs: the Python property tables (cell volume–elevation curves, face
area–elevation profiles, cell minimum elevations) are built by sampling terrain
pixel centres, and RASMapper builds them differently. On these examples cell
minimum elevations differ by up to 18 ft and face conveyance tables by about
20–35 % (median). Closing that gap means matching RASMapper's property-table
algorithm. Until then, use Workflow A for production runs.

The Workflow C mesh itself is close to RASMapper's: from the same seed points
and perimeter it produces the same interior cells and within 0–4 faces of
RASMapper's count on every example (the difference is in how a few perimeter
corners are split). C's results differ from the GUI for the same reasons as B's.

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
the mesh, samples the terrain GeoTIFF for cell volume–elevation curves and face
area–elevation profiles, samples land cover for Manning's *n* and percent
impervious, reads infiltration parameters if `Infiltration.hdf` is present,
assigns perimeter faces to BC lines, and assembles `p01.tmp.hdf` from the
`.p01` text.

### Workflow C — no g01.hdf

**When:** there is no `g01.hdf`. The preprocessor reads the seed points,
perimeter and BC lines from the `.g01` text, clips each seed's Voronoi cell to
the perimeter, and stitches the cells into a mesh that follows RASMapper's
conventions (perimeter faces between face points, one ghost cell per perimeter
face, counter-clockwise cells). It then continues as Workflow B.

## Repository layout

```
ras_preprocess.py          preprocessor (Workflows A / B / C)
hecras_run.py              runs the engines natively or in Docker, checks the result
run_simulation.sh          zip in → preprocess → simulate → zip out
tools/
  make_workflow_input.py   derive a Workflow B or C input from a GUI project
  compare_results.py       compare max WSE between two plan HDFs
tests/test_pipeline.py     pytest: all workflows, all examples, all engines
Examples/                  six HEC-RAS 6.6 projects (see Examples/README.md)
bin/6.6/  bin/7.0.1/       HEC-RAS Linux engines (git LFS)
libs/                      shared libraries, identical for both versions (git LFS)
Dockerfile                 hecras-linux:<version> image
RAS_v.6.6_Linux.pdf        HEC's documentation for the Linux engines
remove_HDF5_Results_Sed.py strip results from a plan HDF (as Workflow A does)
```

## Tests

```bash
make test-fast                         # preprocessing only, ~1 min
make test                              # all workflows, every engine under bin/
pytest --hecras-version 7.0.1          # one engine
pytest --runner docker                 # through the Docker images
```

Workflow A tests must pass. Workflow B and C tests check that the plan runs;
their accuracy tests are marked as expected failures. The run ends with a
table of volume error and WSE difference for every run.

## Adding another engine version

Put the three executables in `bin/<version>/`. If its `libs/` differ from the
shared ones, add them too and point `hecras_run.py` at them. Then
`pytest --hecras-version <version>` and
`docker build --build-arg HECRAS_VERSION=<version> -t hecras-linux:<version> .`.
