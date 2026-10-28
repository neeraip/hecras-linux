"""Pipeline tests: ras_preprocess.py + the HEC-RAS Linux engines.

What is checked, by workflow
  A  The shipped examples (GUI-computed geometry and plan). The Linux engines
     must reproduce the Windows GUI results: max WSE within 0.005 ft
     (0.01 ft for Muncie) and volume error under 0.5 %.
  B  Each example with its geometry tables stripped (tools/make_workflow_input.py),
     plus BEFORE_RUN as shipped. property_tables.py computes the tables, which
     must equal RASMapper's dataset for dataset; the plan must run. Accuracy
     against the GUI is marked as an expected failure until it is measured on
     every example.

The B runs take their boundary conditions from the GUI plan with
--event-conditions because the examples' DSS files are not in the repo.

    pytest -m "not engine"     # preprocessing only (~1 min)
    pytest                     # everything, every installed engine (~10 min)
"""
from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

import h5py
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

import hecras_run  # noqa: E402
from tools.compare_results import compare_max_wse  # noqa: E402
from tools.make_workflow_input import make_input  # noqa: E402

EXAMPLES = REPO / "Examples"
SCRIPT = REPO / "ras_preprocess.py"

# Example folder -> HEC-RAS project name. All are GUI-computed (Workflow A).
GUI_EXAMPLES = {
    "AFTER_RUN": "test_hdf",
    "BEC": "BEC",
    "BEC_WO_Infiltration": "model1",
    "VA": "VA",
    "Muncie": "Muncie",
}
WSE_TOL_FT = {"Muncie": 0.01}       # README: Muncie differs by 0.001 ft
DEFAULT_WSE_TOL_FT = 0.005
VOLUME_ERR_TOL_PCT = 0.5
# The examples finish in about a minute; an unstable B model can spin for hours.
ENGINE_TIMEOUT_S = 900

# B's tables match RASMapper's (test_workflow_b_tables_match_rasmapper), but
# only BEC and BEC_WO_Infiltration have been run through the engine since.
ACCURACY_XFAIL = pytest.mark.xfail(
    reason="engine accuracy not yet re-measured on every example; see README "
           "'Workflow B accuracy'", strict=False)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def preprocess(project_dir: Path, out_dir: Path, project: str,
               event_conditions: Path | None = None) -> str:
    """Run ras_preprocess.py; return its log. Raises on failure."""
    cmd = [sys.executable, str(SCRIPT), str(project_dir), "--output-dir", str(out_dir),
           "--project-name", project]
    if event_conditions:
        cmd += ["--event-conditions", str(event_conditions)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    log = r.stdout + r.stderr
    if r.returncode != 0:
        raise RuntimeError(f"ras_preprocess.py failed:\n{log[-4000:]}")
    return log


def detected_workflow(log: str) -> str:
    m = re.search(r"=== WORKFLOW (\w)", log)
    return m.group(1) if m else "?"


def stage_bc_files(src_dir: Path, out_dir: Path, project: str) -> None:
    """Use the GUI's .b01/.x01 where the example has them."""
    for ext in (".b01", ".x01"):
        src = src_dir / f"{project}{ext}"
        if src.exists():
            shutil.copy2(src, out_dir / src.name)


def area_group(hdf: h5py.File) -> h5py.Group:
    areas = hdf["Geometry/2D Flow Areas"]
    return areas[[k for k in areas if isinstance(areas[k], h5py.Group)][0]]


def has_property_tables(g01_hdf: Path) -> bool:
    with h5py.File(g01_hdf, "r") as f:
        return "Cells Volume Elevation Info" in area_group(f)


_RUNS: dict = {}


@pytest.fixture(scope="session")
def run_dir(tmp_path_factory):
    return tmp_path_factory.mktemp("runs")


def run_case(base: Path, example: str, workflow: str, engine) -> tuple:
    """Prepare and run one example through one workflow, once per session.

    Returns (result, reference plan HDF). A failure is cached and re-raised.
    """
    key = (example, workflow, engine)
    if key not in _RUNS:
        tmp_path = base / f"{example}-{workflow}-{engine[0]}"
        tmp_path.mkdir()
        try:
            _RUNS[key] = _run_case(tmp_path, example, workflow, engine)
        except Exception as e:          # noqa: BLE001 - re-raised for every caller
            _RUNS[key] = e
    if isinstance(_RUNS[key], Exception):
        raise _RUNS[key]
    return _RUNS[key]


def _run_case(tmp_path: Path, example: str, workflow: str, engine) -> tuple:
    version, runner = engine
    project = GUI_EXAMPLES[example]
    src = EXAMPLES / example
    reference = src / f"{project}.p01.hdf"
    out = tmp_path / "run"
    if workflow == "A":
        log = preprocess(src, out, project)
    else:
        inp = tmp_path / "input"
        make_input(src, inp, workflow)
        log = preprocess(inp, out, project, event_conditions=reference)
    assert detected_workflow(log) == workflow
    stage_bc_files(src, out, project)
    result = hecras_run.run_engine(out, project, version=version, runner=runner,
                                   timeout=ENGINE_TIMEOUT_S)
    return result, reference


# ---------------------------------------------------------------------------
# Preprocessing only (no engine)
# ---------------------------------------------------------------------------

class TestPreprocess:

    @pytest.mark.parametrize("example", GUI_EXAMPLES)
    def test_workflow_a_copies_geometry_and_strips_results(self, tmp_path, example):
        project = GUI_EXAMPLES[example]
        log = preprocess(EXAMPLES / example, tmp_path, project)
        assert detected_workflow(log) == "A"
        assert (tmp_path / f"{project}.g01.hdf").exists()
        with h5py.File(tmp_path / f"{project}.p01.tmp.hdf", "r") as f:
            assert "Results" not in f

    def test_before_run_is_workflow_b(self, tmp_path):
        log = preprocess(EXAMPLES / "BEFORE_RUN", tmp_path, "test_hdf")
        assert detected_workflow(log) == "B"
        assert has_property_tables(tmp_path / "test_hdf.g01.hdf")

    @pytest.mark.parametrize("example", GUI_EXAMPLES)
    def test_workflow_b_writes_property_tables(self, tmp_path, example):
        project = GUI_EXAMPLES[example]
        make_input(EXAMPLES / example, tmp_path / "in", "B")
        log = preprocess(tmp_path / "in", tmp_path / "out", project)
        assert detected_workflow(log) == "B"
        g01 = tmp_path / "out" / f"{project}.g01.hdf"
        assert has_property_tables(g01)
        with h5py.File(g01, "r") as f:
            ag = area_group(f)
            n_cells = len(ag["Cells Center Coordinate"])
            n_faces = len(ag["Faces Cell Indexes"])
            assert len(ag["Cells Volume Elevation Info"]) == n_cells
            assert len(ag["Faces Area Elevation Info"]) == n_faces
            assert len(f["Geometry/Boundary Condition Lines/External Faces"]) > 0

    @pytest.mark.parametrize("example", GUI_EXAMPLES)
    def test_workflow_b_tables_match_rasmapper(self, tmp_path, example):
        """On RASMapper's mesh, every dataset RASMapper writes for the 2D area
        must come out identical, byte for byte."""
        project = GUI_EXAMPLES[example]
        make_input(EXAMPLES / example, tmp_path / "in", "B")
        preprocess(tmp_path / "in", tmp_path / "out", project)
        with h5py.File(tmp_path / "out" / f"{project}.g01.hdf", "r") as ours, \
             h5py.File(EXAMPLES / example / f"{project}.g01.hdf", "r") as gui:
            a, g = area_group(ours), area_group(gui)
            names = []
            g.visititems(lambda n, o: names.append(n) if isinstance(o, h5py.Dataset) else None)
            for name in names:
                assert name in a, f"{name} missing"
                assert a[name][()].tobytes() == g[name][()].tobytes(), f"{name} differs"


def f_attr(hdf: h5py.File, field: str):
    return hdf["Geometry/2D Flow Areas/Attributes"][()][field][0]


# ---------------------------------------------------------------------------
# Full pipeline (needs the engines)
# ---------------------------------------------------------------------------

@pytest.mark.engine
class TestWorkflowA:

    @pytest.mark.parametrize("example", GUI_EXAMPLES)
    def test_reproduces_gui_results(self, run_dir, example, engine):
        result, reference = run_case(run_dir, example, "A", engine)
        assert result.volume_error_pct is not None
        assert abs(result.volume_error_pct) < VOLUME_ERR_TOL_PCT
        cmp = compare_max_wse(result.plan_hdf, reference)
        tol = WSE_TOL_FT.get(example, DEFAULT_WSE_TOL_FT)
        assert cmp.max_abs <= tol, f"max WSE diff {cmp.max_abs:.4f} ft > {tol} ft"


@pytest.mark.engine
@pytest.mark.experimental
class TestWorkflowB:

    @pytest.mark.parametrize("example", GUI_EXAMPLES)
    def test_runs(self, run_dir, example, engine):
        result, _ = run_case(run_dir, example, "B", engine)
        assert result.volume_error_pct is not None

    def test_before_run_runs(self, tmp_path, engine):
        version, runner = engine
        out = tmp_path / "run"
        preprocess(EXAMPLES / "BEFORE_RUN", out, "test_hdf",
                   event_conditions=EXAMPLES / "AFTER_RUN" / "test_hdf.p01.hdf")
        stage_bc_files(EXAMPLES / "AFTER_RUN", out, "test_hdf")
        result = hecras_run.run_engine(out, "test_hdf", version=version, runner=runner)
        assert result.volume_error_pct is not None

    @ACCURACY_XFAIL
    @pytest.mark.parametrize("example", GUI_EXAMPLES)
    def test_matches_gui(self, run_dir, example, engine):
        result, reference = run_case(run_dir, example, "B", engine)
        cmp = compare_max_wse(result.plan_hdf, reference)
        assert abs(result.volume_error_pct) < VOLUME_ERR_TOL_PCT
        assert cmp.max_abs <= 0.1, f"max WSE diff {cmp.max_abs:.3f} ft"
