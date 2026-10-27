#!/usr/bin/env python3
"""Derive a Workflow B or C input project from a GUI-computed one.

    python3 tools/make_workflow_input.py Examples/BEC /tmp/BEC_B --workflow B
    python3 tools/make_workflow_input.py Examples/BEC /tmp/BEC_C --workflow C

  B  keeps the RASMapper mesh in <name>.g01.hdf but removes everything that
     "Compute Geometry" adds (hydraulic property tables, BC external faces,
     GeomPreprocess output), as if the mesh had just been drawn.
  C  removes <name>.g01.hdf altogether, leaving the seed points and perimeter
     in the <name>.g01 text file.

Both remove the plan results (<name>.p01.hdf) so ras_preprocess.py cannot
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
    "Infiltration", "Percent Impervious",
]
GEOMETRY_TABLES = [
    "Geometry/Boundary Condition Lines/External Faces",
    "Geometry/Cross Sections",
    "Geometry/GeomPreprocess",
    "Geometry/Structures/Property Tables",
]


def project_name(project_dir: Path) -> str:
    g01 = sorted(project_dir.glob("*.g01"))
    if not g01:
        raise SystemExit(f"No .g01 file in {project_dir}")
    return g01[0].stem


def strip_geometry_tables(g01_hdf: Path) -> None:
    with h5py.File(g01_hdf, "a") as f:
        areas = f["Geometry/2D Flow Areas"]
        for area in [k for k in areas if isinstance(areas[k], h5py.Group)]:
            for name in AREA_TABLES:
                if name in areas[area]:
                    del areas[area][name]
        for path in GEOMETRY_TABLES:
            if path in f:
                del f[path]


def make_input(src: Path, dst: Path, workflow: str, plan: str = "01") -> str:
    workflow = workflow.upper()
    if workflow not in ("B", "C"):
        raise ValueError("workflow must be B or C")
    if dst.exists():
        raise FileExistsError(dst)
    name = project_name(src)
    shutil.copytree(src, dst)
    for f in dst.glob(f"{name}.p{plan}*.hdf"):
        f.unlink()
    g01_hdf = dst / f"{name}.g{plan}.hdf"
    if workflow == "B":
        strip_geometry_tables(g01_hdf)
    elif g01_hdf.exists():
        g01_hdf.unlink()
    return name


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("src", type=Path, help="GUI-computed HEC-RAS project folder")
    ap.add_argument("dst", type=Path, help="new folder to create")
    ap.add_argument("--workflow", required=True, choices=["B", "C", "b", "c"])
    args = ap.parse_args(argv)
    name = make_input(args.src, args.dst, args.workflow)
    print(f"Workflow {args.workflow.upper()} input for {name}: {args.dst}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
