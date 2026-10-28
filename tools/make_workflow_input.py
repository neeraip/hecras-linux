#!/usr/bin/env python3
"""Derive a Workflow B input project from a GUI-computed one.

    python3 tools/make_workflow_input.py Examples/BEC /tmp/BEC_B --workflow B

It keeps the RASMapper mesh in <name>.g01.hdf but removes everything that
"Compute Geometry" adds (hydraulic property tables, BC external faces,
GeomPreprocess output), as if the mesh had just been drawn, and removes the
plan results (<name>.p01.hdf) so ras_preprocess.py cannot
fall back to Workflow A. The GUI results are the reference to compare
against, so keep the original project.
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import h5py

# Datasets/groups RASMapper adds to a 2D flow area when it computes geometry.
AREA_TABLES = [
    "Cells Center Manning's n", "Cells Minimum Elevation", "Cells Surface Area",
    "Cells Volume Elevation Info", "Cells Volume Elevation Values",
    "Faces Area Elevation Info", "Faces Area Elevation Values",
    "Faces Low Elevation Centroid", "Faces Minimum Elevation",
    "Faces Profile Info", "Faces Profile Values",
    "Infiltration", "Percent Impervious",
]
GEOMETRY_TABLES = [
    "Geometry/Boundary Condition Lines/External Faces",
    "Geometry/Cross Sections",
    "Geometry/GeomPreprocess",
    "Geometry/Structures/Property Tables",
]


def project_name(project_dir: Path) -> str:
    """The HEC-RAS project: a .prj with a plan beside it (a folder can hold
    other geometries, such as a Backup.g01, and coordinate-system .prj files)."""
    for prj in sorted(project_dir.glob("*.prj")):
        if any(project_dir.glob(f"{prj.stem}.p[0-9][0-9]")):
            return prj.stem
    g01 = sorted(project_dir.glob("*.g01"))
    if not g01:
        raise SystemExit(f"No HEC-RAS project in {project_dir}")
    return g01[0].stem


def plan_geometry(plan: Path, default: str) -> str:
    """The geometry a plan file names ("Geom File=g08" gives "g08")."""
    for line in plan.read_text(errors="replace").splitlines():
        if line.startswith("Geom File=") and line.split("=", 1)[1].strip():
            return line.split("=", 1)[1].strip()
    return default


def strip_geometry_tables(g01_hdf: Path) -> None:
    with h5py.File(g01_hdf, "a") as f:
        if "Geometry/2D Flow Areas" not in f:
            raise SystemExit(f"{g01_hdf.name} has no 2D flow area")
        areas = f["Geometry/2D Flow Areas"]
        for area in [k for k in areas if isinstance(areas[k], h5py.Group)]:
            for name in AREA_TABLES:
                if name in areas[area]:
                    del areas[area][name]
        for path in GEOMETRY_TABLES:
            if path in f:
                del f[path]


def make_input(src: Path, dst: Path, workflow: str, plan: str = "01", name: str = None) -> str:
    workflow = workflow.upper()
    if workflow != "B":
        raise ValueError("workflow must be B")
    if dst.exists():
        raise FileExistsError(dst)
    name = name or project_name(src)
    shutil.copytree(src, dst)
    for f in dst.glob(f"{name}.p{plan}*.hdf"):
        f.unlink()
    geom = plan_geometry(dst / f"{name}.p{plan}", f"g{plan}")
    strip_geometry_tables(dst / f"{name}.{geom}.hdf")
    return name


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("src", type=Path, help="GUI-computed HEC-RAS project folder")
    ap.add_argument("dst", type=Path, help="new folder to create")
    ap.add_argument("--workflow", required=True, choices=["B", "b"])
    ap.add_argument("--plan", default="01", help="plan number (default 01)")
    ap.add_argument("--project", help="project name, when the folder holds more than one")
    args = ap.parse_args(argv)
    name = make_input(args.src, args.dst, args.workflow, args.plan.zfill(2), args.project)
    print(f"Workflow {args.workflow.upper()} input for {name}: {args.dst}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
