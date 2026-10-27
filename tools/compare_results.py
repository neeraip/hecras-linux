#!/usr/bin/env python3
"""Compare maximum water surface elevations between two HEC-RAS plan HDFs.

    python3 tools/compare_results.py result.p01.hdf reference.p01.hdf

Cells are matched by index when both runs share a mesh, otherwise each
reference cell is matched to the nearest result cell centre (for Workflow C,
whose mesh is rebuilt). All active reference cells are compared, dry ones
included, so a run that floods a cell the reference leaves dry still counts.
"""
from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np

SUMMARY = "Results/Unsteady/Output/Output Blocks/Base Output/Summary Output/2D Flow Areas"


@dataclass
class WseComparison:
    max_abs: float
    p99_abs: float
    mean_abs: float
    n_cells: int
    same_mesh: bool


def _area(f: h5py.File) -> str:
    return [k for k in f[SUMMARY] if isinstance(f[SUMMARY][k], h5py.Group)][0]


def _max_wse(f: h5py.File, area: str) -> np.ndarray:
    return f[f"{SUMMARY}/{area}/Maximum Water Surface"][0].astype(np.float64)


def _centres(f: h5py.File, area: str) -> np.ndarray:
    return f[f"Geometry/2D Flow Areas/{area}/Cells Center Coordinate"][()]


def compare_max_wse(result_hdf: Path, reference_hdf: Path) -> WseComparison:
    with h5py.File(result_hdf, "r") as rf, h5py.File(reference_hdf, "r") as ff:
        ra, fa = _area(rf), _area(ff)
        r, ref = _max_wse(rf, ra), _max_wse(ff, fa)
        rc, fc = _centres(rf, ra)[: len(r)], _centres(ff, fa)[: len(ref)]
    same = rc.shape == fc.shape and np.allclose(rc, fc)
    if not same:
        from scipy.spatial import cKDTree
        r = r[cKDTree(rc).query(fc)[1]]
    active = np.isfinite(ref) & (ref < 1e10) & np.isfinite(r)
    if not active.any():
        return WseComparison(0.0, 0.0, 0.0, 0, same)
    d = np.abs(r[active] - ref[active])
    return WseComparison(float(d.max()), float(np.percentile(d, 99)), float(d.mean()),
                         int(active.sum()), same)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("result", type=Path)
    ap.add_argument("reference", type=Path)
    a = ap.parse_args(argv)
    c = compare_max_wse(a.result, a.reference)
    how = "same mesh" if c.same_mesh else "nearest cell"
    print(f"{c.n_cells} cells ({how}): max |dWSE| {c.max_abs:.4f}, "
          f"99th pct {c.p99_abs:.4f}, mean {c.mean_abs:.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
