"""2D flow area property tables, computed the way RASMapper computes them.

A port of RASMapper's "Compute Property Tables" (RasMapperLib, HEC-RAS
6.6/7.0): face profiles and area-elevation tables, cell volume-elevation
tables, cell minimum elevations and face low-elevation centroids. Given the
same mesh, terrain and settings it reproduces the tables the GUI writes.

The arithmetic follows RASMapper's, including where it rounds to 32-bit
floats, because the tables are compared value for value with the GUI's.

Terrain is a RASMapper terrain GeoTIFF (tiled, with its own pyramid levels).
Values are interpolated in "render blocks": the square between four pixel
centres, split into four triangles that meet at the block centre, whose value
is the mean of the four corners.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
from osgeo import gdal

from . import tin

F32 = np.float32
NODATA = -9999.0
# .NET's float.NaN, which RASMapper writes: a NaN with the sign bit set.
NAN = np.frombuffer(b"\x00\x00\xc0\xff", dtype=F32)[0]

# RASMapper's default settings (RASD2FlowArea), used when a geometry has none.
DEFAULTS = {
    "cell_volume_tolerance": 0.01,
    "cell_min_area_fraction": 0.01,
    "face_profile_tolerance": 0.01,
    "face_area_elevation_tolerance": 0.01,
    "face_conveyance_ratio": 0.02,
    "laminar_depth": 0.2,
}

CELL_RASTER_SIZE = 10000      # pixels in the raster each cell is resampled onto
HISTOGRAM_BINS = 256
CONVEYANCE_FACTOR = F32(1.486)  # US customary units


# ----------------------------------------------------------------------------
# Terrain
# ----------------------------------------------------------------------------

@dataclass
class _Level:
    z: np.ndarray        # float32 pixel values, NaN for no data
    cell: float          # pixel size
    x0: float            # x of the first pixel centre (render-block origin)
    y0: float            # y of the first pixel centre
    tile: int            # TIFF tile size in pixels
    mask: np.ndarray = None  # render-block flags from the terrain HDF, if any:
                             # 1 has data; 2, 4, 8, 16 its bottom, right, left,
                             # top edge is a stitch (borders no data)
    rb_x0: float = 0.0       # RasterFileInfo.RenderBlockTiles: the grid of
    rb_y0: float = 0.0       # tiles of render blocks, whose spacing RASMapper
    rb_cell: float = 0.0     # recomputes from its extent
    ntw: int = 0             # tiles wide and tall
    ntt: int = 0

    def tile_has_data(self, tr, tc):
        """RasterFileInfo.TileHasAnyData."""
        t = self.tile
        if not (0 <= tr < self.ntt and 0 <= tc < self.ntw):
            return False
        if self.mask is not None:
            return bool((self.mask[tr * t:(tr + 1) * t, tc * t:(tc + 1) * t] & 1).any())
        return bool(np.isfinite(self.z[tr * t:(tr + 1) * t, tc * t:(tc + 1) * t]).any())

    def tile_frame(self, tr, tc):
        """RasterM.SetCellTileRasterM: a tile's left x, top y and block size,
        each recomputed from the tile's own corners as RASMapper does."""
        min_x = self.rb_x0 + self.rb_cell * tc
        max_x = self.rb_x0 + self.rb_cell * (tc + 1)
        max_y = self.rb_y0 - self.rb_cell * tr
        min_y = self.rb_y0 - self.rb_cell * (tr + 1)
        cell = np.maximum((max_x - min_x) / self.tile, (max_y - min_y) / self.tile)
        return min_x, max_y, cell

    def block_bounds(self, r, c):
        """Render block (r, c) as min x, max x, max y, min y, centre x, centre y."""
        t = self.tile
        min_x, max_y, cell = self.tile_frame(r // t, c // t)
        lr, lc = r % t, c % t
        return (min_x + cell * lc, min_x + cell * (lc + 1), max_y - cell * lr,
                max_y - cell * (lr + 1), min_x + cell * (lc + 0.5), max_y - cell * (lr + 0.5))

    def locate(self, x, y, tr=None, tc=None):
        """The render block holding (x, y) and the fractions across it
        (TerrainLayer.ResampleInternal / RasterM.ColAndFraction, RowAndFraction),
        in the given tile or the one the point falls in."""
        t = self.tile
        if tc is None:
            tc = np.floor((x - self.rb_x0) / self.rb_cell).astype(np.int64)
        if tr is None:
            tr = np.floor((self.rb_y0 - y) / self.rb_cell).astype(np.int64)
        min_x, max_y, cell = self.tile_frame(tr, tc)
        lc = np.floor((x - min_x) / cell).astype(np.int64)
        lr = np.floor((max_y - y) / cell).astype(np.int64)
        xf = ((x - (min_x + cell * lc)) / cell).astype(F32)
        yf = (((max_y - cell * lr) - y) / cell).astype(F32)
        xf = np.where(lc < 0, F32(0), np.where(lc >= t, F32(1), xf)).astype(F32)
        yf = np.where(lr < 0, F32(0), np.where(lr >= t, F32(1), yf)).astype(F32)
        lc, lr = np.clip(lc, 0, t - 1), np.clip(lr, 0, t - 1)
        return tr * t + lr, tc * t + lc, xf, yf


class Terrain:
    """A RASMapper terrain raster and its pyramid levels."""

    def __init__(self, path: str):
        ds = gdal.Open(path)
        band = ds.GetRasterBand(1)
        gt = ds.GetGeoTransform()
        nodata = band.GetNoDataValue()
        tile = band.GetBlockSize()[0]
        ext_min_x, ext_max_y = gt[0], gt[3]
        ext_max_x = gt[0] + ds.RasterXSize * gt[1]
        ext_min_y = gt[3] + ds.RasterYSize * gt[5]
        cell = max((ext_max_x - ext_min_x) / ds.RasterXSize,
                   (ext_max_y - ext_min_y) / ds.RasterYSize)
        bands = [band] + [band.GetOverview(i) for i in range(band.GetOverviewCount())]
        self.levels = []
        for i, b in enumerate(bands):
            z = b.ReadAsArray().astype(F32)
            if nodata is not None:
                z[z == F32(nodata)] = np.nan
            if i:
                cell *= 2.0
            t = b.GetBlockSize()[0] or tile
            ntw, ntt = -(-z.shape[1] // t), -(-z.shape[0] // t)
            span = t * cell
            x, y = ext_min_x + cell * 0.5, ext_max_y - cell * 0.5
            rb_cell = max(((x + ntw * span) - x) / ntw, (y - (y - ntt * span)) / ntt)
            self.levels.append(_Level(z, cell, x, y, t, None, x, y, rb_cell, ntw, ntt))
        self.min_x, self.max_y = gt[0], gt[3]
        self.max_x = gt[0] + ds.RasterXSize * gt[1]
        self.min_y = gt[3] + ds.RasterYSize * gt[5]
        ds = None
        self._load_masks(Path(path))

    def _load_masks(self, tif: Path):
        """The render-block masks RASMapper keeps in the terrain HDF next to the
        raster (Terrain/<layer>/<level>/Mask, one 256x256 tile per row)."""
        for hdf in tif.parent.glob("*.hdf"):
            try:
                hf = h5py.File(str(hdf), "r")
            except OSError:
                continue
            with hf:
                if "Terrain" not in hf:
                    continue
                for name, grp in hf["Terrain"].items():
                    if not isinstance(grp, h5py.Group) or grp.attrs.get("File", b"").decode() != tif.name:
                        continue
                    for i, lv in enumerate(self.levels):
                        if f"{i}/Mask" not in grp:
                            continue
                        tiles = grp[f"{i}/Mask"][()]
                        t = lv.tile
                        nr, nc = lv.z.shape
                        ntw = -(-nc // t)
                        ntt = len(tiles) // ntw
                        m = tiles.reshape(ntt, ntw, t, t).transpose(0, 2, 1, 3).reshape(ntt * t, ntw * t)
                        lv.mask = m[:nr, :nc]
                    return

    def level_for(self, cell_size: float) -> int:
        """RasterFileInfo.ComputeLevel: the first level at least as coarse."""
        for i, lv in enumerate(self.levels):
            if abs(lv.cell - cell_size) < 1e-5 or lv.cell > cell_size:
                return i
        return len(self.levels) - 1

    @staticmethod
    def _corners(lv: _Level, r, c):
        """Corner values a (top-left), b, c, d of render blocks (r, c); NaN if
        the block is off the raster or touches no data."""
        nr, nc = lv.z.shape
        ok = (r >= 0) & (r < nr - 1) & (c >= 0) & (c < nc - 1)
        rr, cc = np.where(ok, r, 0), np.where(ok, c, 0)
        a = lv.z[rr, cc]
        b = lv.z[rr, cc + 1]
        cv = lv.z[rr + 1, cc + 1]
        d = lv.z[rr + 1, cc]
        if lv.mask is not None:
            ok &= (lv.mask[rr, cc] & 1) != 0
        ok &= np.isfinite(a) & np.isfinite(b) & np.isfinite(cv) & np.isfinite(d)
        return ok, a, b, cv, d

    def sample(self, x: np.ndarray, y: np.ndarray, level: int = 0, tile=None) -> np.ndarray:
        """Terrain at points (RenderBlock.Interpolate); NaN where there is none."""
        lv = self.levels[level]
        tr, tc = tile if tile is not None else (None, None)
        r, c, xf, yf = lv.locate(np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64), tr, tc)
        ok, a, b, cv, d = self._corners(lv, r, c)
        return np.where(ok, _interpolate(a, b, cv, d, xf, yf), np.nan)


class TerrainLayer:
    """A RASMapper terrain: its terrain HDF's rasters in priority order
    (TerrainLayer.RasterFilesInfo), or a single GeoTIFF."""

    def __init__(self, path: str):
        path = Path(path)
        self.rasters = []
        if path.suffix.lower() != ".hdf":
            self.rasters.append(Terrain(str(path)))
            return
        with h5py.File(str(path), "r") as hf:
            layers = []
            for name, grp in hf["Terrain"].items():
                if isinstance(grp, h5py.Group) and "File" in grp.attrs:
                    layers.append((int(grp.attrs.get("Priority", 0)), grp.attrs["File"].decode()))
        for _, tif in sorted(layers, key=lambda pr: pr[0]):
            self.rasters.append(Terrain(str(path.parent / tif)))
        if not self.rasters:
            raise ValueError(f"{path}: no rasters listed under /Terrain")
        # TerrainLayer.Level0StitchTIN: the triangles RASMapper stores to bridge
        # the seams between rasters, sampled along faces after the rasters.
        with h5py.File(str(path), "r") as hf:
            pts = hf.get("Terrain/Stitch TIN Points")
            tri = hf.get("Terrain/Stitch TIN Triangles")
            if pts is not None and tri is not None and len(pts) and len(tri):
                self.stitch_points = pts[()].astype(np.float64)
                self.stitch_triangles = tri[()].astype(np.int64)
                v = self.stitch_points[self.stitch_triangles][:, :, :2]
                self.stitch_bounds = np.concatenate([v.min(axis=1), v.max(axis=1)], axis=1)

    stitch_points = None
    stitch_triangles = None

    def stitch_profile(self, pts: np.ndarray) -> list:
        """Tin.ComputeProfiles on the stitch TIN for one polyline: segment ends
        inside a triangle and crossings of triangle edges, as (station, z),
        sorted and filtered at zero tolerance."""
        if self.stitch_triangles is None:
            return []
        P, T = self.stitch_points, self.stitch_triangles
        out = []
        start = 0.0
        for (ax, ay), (bx, by) in zip(pts[:-1], pts[1:]):
            length = float(np.hypot(bx - ax, by - ay))
            lo_x, hi_x, lo_y, hi_y = min(ax, bx), max(ax, bx), min(ay, by), max(ay, by)
            b = self.stitch_bounds
            near = np.nonzero((b[:, 0] <= hi_x) & (b[:, 2] >= lo_x) & (b[:, 1] <= hi_y) & (b[:, 3] >= lo_y))[0]
            for t in near:
                p0, p1, p2 = P[T[t, 0]], P[T[t, 1]], P[T[t, 2]]
                for (px, py, s) in ((ax, ay, start), (bx, by, start + length)):
                    n1, n2, n3 = p1[1] - p2[1], p0[0] - p2[0], p2[0] - p1[0]
                    n4 = p0[1] - p2[1]
                    den = n1 * n2 + n3 * n4
                    dx, dy = px - p2[0], py - p2[1]
                    l0 = (n1 * dx + n3 * dy) / den
                    l1 = (n2 * dy - n4 * dx) / den
                    l2 = 1.0 - l0 - l1
                    if 0 <= l0 <= 1 and 0 <= l1 <= 1 and 0 <= l2 <= 1:
                        out.append((F32(s), F32(p0[2] * l0 + p1[2] * l1 + p2[2] * l2)))
                for e0, e1 in ((p0, p1), (p1, p2), (p2, p0)):
                    hit, t1, t2 = _segment_intersections(ax, ay, bx, by, np.array([e0[0]]), np.array([e0[1]]),
                                                         np.array([e1[0]]), np.array([e1[1]]))
                    if hit[0]:
                        dx, dy = bx - ax, by - ay
                        ex, ey = e1[0] - e0[0], e1[1] - e0[1]
                        den = ey * dx - ex * dy
                        u = (ex * (ay - e0[1]) - ey * (ax - e0[0])) / den   # unclamped t1
                        x, y = ax + u * dx, ay + u * dy
                        st = start + float(np.hypot(x - ax, y - ay))
                        out.append((F32(st), F32(e0[2] + (e1[2] - e0[2]) * t2[0])))
            start += length
        if not out:
            return []
        _dotnet_sort(out, key=lambda q: q[0])
        keep = _dp_keep(np.array([q[0] for q in out], F32), np.array([q[1] for q in out], F32), 0.0)
        return [out[i] for i in keep]


def _as_layer(terrain) -> "TerrainLayer":
    if isinstance(terrain, TerrainLayer):
        return terrain
    layer = TerrainLayer.__new__(TerrainLayer)
    layer.rasters = [terrain]
    layer.stitch_triangles = None
    return layer


def _stitch_segments(lv: _Level, r0, r1, c0, c1) -> list:
    """RenderBlock.SetRowCol's stitches for render blocks rows r0..r1, cols
    c0..c1: the flagged edges of blocks with data, as ((x, y, z, m), (x, y, z, m))
    with m the block size, running clockwise round the block (so no data lies
    on their left)."""
    if lv.mask is None:
        return []
    nr, nc = lv.z.shape
    r0, c0 = max(r0, 0), max(c0, 0)
    r1, c1 = min(r1, nr - 2), min(c1, nc - 2)
    if r1 < r0 or c1 < c0:
        return []
    m = lv.mask[r0:r1 + 1, c0:c1 + 1].astype(np.int64)
    rows, cols = np.nonzero(((m & 1) != 0) & ((m & 0x1E) != 0))
    segs = []
    for r, c in zip(rows + r0, cols + c0):
        b = int(lv.mask[r, c])
        minx, maxx, maxy, miny, _, _ = (float(v) for v in lv.block_bounds(r, c))
        m = float(lv.tile_frame(r // lv.tile, c // lv.tile)[2])
        A = (minx, maxy, float(lv.z[r, c]), m)
        B = (maxx, maxy, float(lv.z[r, c + 1]), m)
        C = (maxx, miny, float(lv.z[r + 1, c + 1]), m)
        D = (minx, miny, float(lv.z[r + 1, c]), m)
        if b & 8:
            segs.append((D, A))
        if b & 0x10:
            segs.append((A, B))
        if b & 4:
            segs.append((B, C))
        if b & 2:
            segs.append((C, D))
    return segs


def _interpolate(a, b, c, d, xf, yf):
    """RenderBlock.Interpolate, in float32 and in RASMapper's operation order."""
    with np.errstate(invalid="ignore"):
        e2 = (((a + b) + c) + d) * F32(0.5)
        ny = F32(1) - yf
        top = (a + (b - a) * xf) + (e2 - (a + b)) * yf
        right = (b + (e2 - (b + c)) * (F32(1) - xf)) + (c - b) * yf
        bottom = (d + (c - d) * xf) + (e2 - (c + d)) * ny
        left = (a + (e2 - (d + a)) * xf) + (d - a) * yf
    return np.where(xf > yf, np.where(xf < ny, top, right),
                    np.where(xf > ny, bottom, left)).astype(F32)


# ----------------------------------------------------------------------------
# Douglas-Peucker on (x, y) float32 points (Point2DCollection.FilterKeepList)
# ----------------------------------------------------------------------------

def _dp_keep(px: np.ndarray, py: np.ndarray, tol: float) -> list:
    n = len(px)
    keep = [0] if n > 0 else []
    if n > 1:
        keep.append(n - 1)
    tol = F32(tol)
    stack = [(0, n - 1)]
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        while stack:
            first, last = stack.pop()
            if last - first <= 1:
                continue
            fx, fy, lx, ly = px[first], py[first], px[last], py[last]
            dist = F32(np.sqrt(float(lx - fx) ** 2 + float(ly - fy) ** 2))
            ix, iy = px[first + 1:last], py[first + 1:last]
            num = np.abs(fx * ly + lx * iy + ix * fy - lx * fy - ix * ly - fx * iy) / dist
            # First index of the largest distance, and only if it beats zero.
            best = F32(0)
            idx = -1
            if len(num):
                m = np.nanmax(num) if np.isfinite(num).any() else np.nan
                if m > best:
                    idx = first + 1 + int(np.argmax(num == m))
                    best = m
            if idx >= 0 and best > tol:
                keep.append(idx)
                stack.append((idx, last))
                stack.append((first, idx))
    return sorted(keep)


# ----------------------------------------------------------------------------
# Faces
# ----------------------------------------------------------------------------

def _segment_intersections(ax, ay, bx, by, sx0, sy0, sx1, sy1, eps=1e-6):
    """SegmentM.Intersects of face segment A-B with many segments s0-s1:
    returns (hit mask, t along A-B, t along s)."""
    dx, dy = bx - ax, by - ay
    ex, ey = sx1 - sx0, sy1 - sy0
    # MBROverlaps with epsilon
    mbr = ((np.minimum(sx0, sx1) <= max(ax, bx) + eps) & (np.maximum(sx0, sx1) >= min(ax, bx) - eps)
           & (np.minimum(sy0, sy1) <= max(ay, by) + eps) & (np.maximum(sy0, sy1) >= min(ay, by) - eps))
    den = ey * dx - ex * dy
    with np.errstate(divide="ignore", invalid="ignore"):
        n6 = ay - sy0
        n7 = ax - sx0
        t1 = (ex * n6 - ey * n7) / den
        t2 = (dx * n6 - dy * n7) / den
    hit = mbr & (den != 0) & (t1 > -eps) & (t1 < 1 + eps) & (t2 > -eps) & (t2 < 1 + eps)
    return hit, np.clip(t1, 0, 1), np.clip(t2, 0, 1)


def _dotnet_sort(items: list, key) -> None:
    """List<T>.Sort(Comparison) on .NET Framework 4.7.2 (introsort: insertion
    sort up to 16 items, median-of-three quicksort, heapsort past the depth
    limit), in place. Matches RASMapper's order of equal keys, which the
    profile filter is sensitive to."""
    n = len(items)
    if n < 2:
        return
    k = [key(v) for v in items]

    def cmp(a, b):
        return -1 if a < b else (0 if a == b else 1)

    def swap(i, j):
        items[i], items[j] = items[j], items[i]
        k[i], k[j] = k[j], k[i]

    def swap_if_greater(i, j):
        if i != j and cmp(k[i], k[j]) > 0:
            swap(i, j)

    def insertion(lo, hi):
        for i in range(lo, hi):
            j = i
            tv, tk = items[i + 1], k[i + 1]
            while j >= lo and cmp(tk, k[j]) < 0:
                items[j + 1], k[j + 1] = items[j], k[j]
                j -= 1
            items[j + 1], k[j + 1] = tv, tk

    def down_heap(i, m, lo):
        dv, dk = items[lo + i - 1], k[lo + i - 1]
        while i <= m // 2:
            child = 2 * i
            if child < m and cmp(k[lo + child - 1], k[lo + child]) < 0:
                child += 1
            if not cmp(dk, k[lo + child - 1]) < 0:
                break
            items[lo + i - 1], k[lo + i - 1] = items[lo + child - 1], k[lo + child - 1]
            i = child
        items[lo + i - 1], k[lo + i - 1] = dv, dk

    def heapsort(lo, hi):
        m = hi - lo + 1
        for i in range(m // 2, 0, -1):
            down_heap(i, m, lo)
        for i in range(m, 1, -1):
            swap(lo, lo + i - 1)
            down_heap(1, i - 1, lo)

    def partition(lo, hi):
        mid = lo + (hi - lo) // 2
        swap_if_greater(lo, mid)
        swap_if_greater(lo, hi)
        swap_if_greater(mid, hi)
        pivot = k[mid]
        swap(mid, hi - 1)
        left, right = lo, hi - 1
        while left < right:
            left += 1
            while cmp(k[left], pivot) < 0:
                left += 1
            right -= 1
            while cmp(pivot, k[right]) < 0:
                right -= 1
            if left >= right:
                break
            swap(left, right)
        swap(left, hi - 1)
        return left

    def intro(lo, hi, depth):
        while hi > lo:
            size = hi - lo + 1
            if size <= 16:
                if size == 1:
                    return
                if size == 2:
                    swap_if_greater(lo, hi)
                    return
                if size == 3:
                    swap_if_greater(lo, hi - 1)
                    swap_if_greater(lo, hi)
                    swap_if_greater(hi - 1, hi)
                    return
                insertion(lo, hi)
                return
            if depth == 0:
                heapsort(lo, hi)
                return
            depth -= 1
            p = partition(lo, hi)
            intro(p + 1, hi, depth)
            hi = p - 1

    intro(0, n - 1, 2 * (n.bit_length() - 1))


def _clip(ax, ay, bx, by, min_x, max_x, min_y, max_y):
    """SegmentM.ClipToExtent: the segment clipped to an extent, or None."""
    if (max(ax, bx) < min_x or min(ax, bx) > max_x or max(ay, by) < min_y or min(ay, by) > max_y):
        return None
    t0, t1 = 0.0, 1.0
    dx, dy = bx - ax, by - ay
    for p, q in ((-dx, ax - min_x), (dx, max_x - ax), (-dy, ay - min_y), (dy, max_y - ay)):
        if p < 0:
            r = q / p
            if r > t1:
                return None
            if r > t0:
                t0 = r
        elif p > 0:
            r = q / p
            if r < t0:
                return None
            if r < t1:
                t1 = r
        elif q < 0:
            return None
    nbx, nby = (ax + t1 * dx, ay + t1 * dy) if t1 < 1.0 else (bx, by)
    nax, nay = (ax + t0 * dx, ay + t0 * dy) if t0 > 0.0 else (ax, ay)
    return nax, nay, nbx, nby


def _walk_cells(ax, ay, bx, by, min_x, max_x, min_y, max_y, cell, rows, cols, tol=0.0):
    """RasterizeSegment.ComputeCells: the raster cells a segment passes
    through, in order from A, as (row, col)."""
    inside = lambda x, y: (min_x - tol <= x <= max_x + tol) and (min_y - tol <= y <= max_y + tol)
    if not (inside(ax, ay) and inside(bx, by)):
        c = _clip(ax, ay, bx, by, min_x, max_x, min_y, max_y)
        if c is None:
            return []
        ax, ay, bx, by = c
    row = lambda y: int(np.floor((max_y - y) / cell))
    col = lambda x: int(np.floor((x - min_x) / cell))
    r0, c0, r1, c1 = row(ay), col(ax), row(by), col(bx)
    r0 = 0 if r0 < 0 else (rows - 1 if r0 == rows else r0)
    c0 = 0 if c0 < 0 else (cols - 1 if c0 == cols else c0)
    r1 = 0 if r1 < 0 else (rows - 1 if r1 == rows else r1)
    c1 = 0 if c1 < 0 else (cols - 1 if c1 == cols else c1)
    dc, dr = c1 - c0, r1 - r0
    dx, dy = bx - ax, by - ay
    if dc == 0 and dr == 0:
        return [(r0, c0)]
    if dc == 0:
        s = 1 if dr > 0 else -1
        return [(r, c0) for r in range(r0, r1 + s, s)]
    if dr == 0:
        s = 1 if dc > 0 else -1
        return [(r0, c) for c in range(c0, c1 + s, s)]
    out = [(r0, c0)]
    r, c = r0, c0
    min_y_of = lambda r: max_y - cell * (r + 1)
    max_y_of = lambda r: max_y - cell * r
    max_x_of = lambda c: min_x + cell * (c + 1)
    min_x_of = lambda c: min_x + cell * c
    while not (c == c1 and r == r1):
        ty = ((min_y_of(r) if dr > 0 else max_y_of(r)) - ay) / dy
        tx = ((max_x_of(c) if dc > 0 else min_x_of(c)) - ax) / dx
        step_row = (ty < tx) and (r < r1 if dr > 0 else r > r1)
        if step_row or c == c1:
            r += 1 if dr > 0 else -1
        else:
            c += 1 if dc > 0 else -1
        out.append((r, c))
    return out


def _tile_profile(terrain, lv, tr, tc, ax, ay, bx, by, start, length):
    """ComputeProfilesOnTile: the segment's points on one tile, in RASMapper's
    order: ends that lie on the tile, then each block it crosses, eight edges
    each (centre to corners, then the sides)."""
    t = lv.tile
    min_x, max_y, cell = lv.tile_frame(tr, tc)
    max_x = lv.rb_x0 + lv.rb_cell * (tc + 1)
    min_y = lv.rb_y0 - lv.rb_cell * (tr + 1)
    pts = []
    for (px, py, s) in ((ax, ay, start), (bx, by, start + length)):
        if min_x - 1e-7 <= px <= max_x + 1e-7 and min_y - 1e-7 <= py <= max_y + 1e-7:
            z = terrain.sample(np.array([px]), np.array([py]), 0, (tr, tc))[0]
            if np.isfinite(z):
                pts.append((F32(s), F32(z)))
    cells = _walk_cells(ax, ay, bx, by, min_x, max_x, min_y, max_y, cell, t, t, 1e-7)
    if not cells:
        return pts
    lr = np.array([rc[0] for rc in cells])
    lc = np.array([rc[1] for rc in cells])
    ok, a, b, c, d = Terrain._corners(lv, tr * t + lr, tc * t + lc)
    lr, lc, a, b, c, d = lr[ok], lc[ok], a[ok], b[ok], c[ok], d[ok]
    if len(lr) == 0:
        return pts
    bminx, bmaxx = min_x + cell * lc, min_x + cell * (lc + 1)
    bmaxy, bminy = max_y - cell * lr, max_y - cell * (lr + 1)
    cx, cy = min_x + cell * (lc + 0.5), max_y - cell * (lr + 0.5)
    ez = ((((a + b) + c) + d) * F32(0.5)).astype(np.float64) * 0.5
    A = (bminx, bmaxy, a.astype(np.float64))
    B = (bmaxx, bmaxy, b.astype(np.float64))
    C = (bmaxx, bminy, c.astype(np.float64))
    D = (bminx, bminy, d.astype(np.float64))
    E = (cx, cy, ez)
    edges = [(E, A), (E, B), (E, C), (E, D), (A, B), (B, C), (C, D), (D, A)]
    hits, ss, zz = [], [], []
    for p, q in edges:
        h, t1, t2 = _segment_intersections(ax, ay, bx, by, p[0], p[1], q[0], q[1])
        hits.append(h)
        ss.append((start + length * t1).astype(F32))
        zz.append((p[2] + (q[2] - p[2]) * t2).astype(F32))
    hits, ss, zz = np.stack(hits, 1), np.stack(ss, 1), np.stack(zz, 1)   # block x edge
    for i in range(len(lr)):
        for e in range(8):
            if hits[i, e]:
                pts.append((ss[i, e], zz[i, e]))
    return pts


def face_profile(terrain, pts: np.ndarray, tol: float) -> tuple:
    """TerrainLayer.ComputeProfiles for one face polyline: (station, elevation)
    as float32 arrays, thinned with the face profile tolerance. Each raster
    adds the points where it has data; the whole is then sorted and thinned."""
    layer = _as_layer(terrain)
    profile = []
    for raster in layer.rasters:
        profile.extend(_raster_profile(raster, pts, tol))
    stitched = layer.stitch_profile(pts) if isinstance(layer, TerrainLayer) else []
    if stitched:
        profile.extend(stitched)
        _dotnet_sort(profile, key=lambda p: p[0])
    if not profile:
        return np.empty(0, F32), np.empty(0, F32)
    total = float(np.sum(np.hypot(np.diff(pts[:, 0]), np.diff(pts[:, 1]))))
    _dotnet_sort(profile, key=lambda p: p[0])
    if abs(float(profile[-1][0]) - total) > 0.001:
        profile.append((F32(total), profile[-1][1]))
    if abs(float(profile[0][0])) > 0.001:
        profile.append((F32(0), profile[0][1]))
    _dotnet_sort(profile, key=lambda p: p[0])
    s = np.array([p[0] for p in profile], F32)
    z = np.array([p[1] for p in profile], F32)
    if tol > 0:
        keep = _dp_keep(s, z, tol)
        s, z = s[keep], z[keep]
    return s, z


def _raster_profile(terrain: Terrain, pts: np.ndarray, tol: float) -> list:
    """TerrainLayer.ComputeProfilesOnRaster: one raster's points along a face
    polyline, each tile's thinned on its own."""
    lv = terrain.levels[0]
    tiles_max_x = lv.rb_x0 + lv.rb_cell * lv.ntw
    tiles_min_y = lv.rb_y0 - lv.rb_cell * lv.ntt
    # Which tiles each segment crosses, in the order RASMapper finds them.
    work = {}
    start = 0.0
    for (ax, ay), (bx, by) in zip(pts[:-1], pts[1:]):
        length = float(np.hypot(bx - ax, by - ay))
        clipped = _clip(ax, ay, bx, by, lv.rb_x0, tiles_max_x, tiles_min_y, lv.rb_y0)
        if clipped is not None:
            for tr, tc in _walk_cells(*clipped, lv.rb_x0, tiles_max_x, tiles_min_y, lv.rb_y0,
                                      lv.rb_cell, lv.ntt, lv.ntw):
                if lv.tile_has_data(tr, tc):
                    work.setdefault((tr, tc), []).append((ax, ay, bx, by, start, length))
        start += length
    profile = []
    for (tr, tc), segs in work.items():
        for ax, ay, bx, by, s0, length in segs:
            part = _tile_profile(terrain, lv, tr, tc, ax, ay, bx, by, s0, length)
            _dotnet_sort(part, key=lambda p: p[0])
            if tol > 0 and part:
                keep = _dp_keep(np.array([p[0] for p in part], F32),
                                np.array([p[1] for p in part], F32), tol)
                part = [part[i] for i in keep]
            profile.extend(part)
    return profile


class _HP:
    """ProfileHydraulicPropertiesTable.HydraulicProperties (float32 fields)."""
    __slots__ = ("elev", "area", "tw", "wp", "cx", "n", "gri")

    def __init__(self, elev, area=0.0, tw=0.0, wp=0.0, cx=0.0, n=0.0, gri=0):
        self.elev, self.area, self.tw, self.wp = F32(elev), F32(area), F32(tw), F32(wp)
        self.cx, self.n, self.gri = F32(cx), F32(n), gri

    def copy(self):
        return _HP(self.elev, self.area, self.tw, self.wp, self.cx, self.n, self.gri)

    def hydraulic_radius(self):
        return F32(0) if self.wp == 0 else F32(self.area / self.wp)

    def conveyance(self):
        n = self.n if self.n > 0 else F32(0.03)
        return F32(float(CONVEYANCE_FACTOR * self.area)
                   * float(np.power(np.float64(self.hydraulic_radius()), 0.6666666865348816)) / float(n))

    def projected_conveyance(self, elev_above):
        n = self.n if self.n > 0 else F32(0.03)
        dz = F32(F32(elev_above) - self.elev)
        area = F32(self.area + F32(self.tw * dz))
        wp = F32(self.wp + dz * F32(2)) if self.gri < 2 else F32(self.wp + dz * F32(self.gri))
        hr = F32(area / wp)
        return F32(float(CONVEYANCE_FACTOR * area)
                   * float(np.power(np.float64(hr), 0.6666666865348816)) / float(n))

    def midpoint(self, hp):
        r = _HP(0)
        r.elev = F32((float(self.elev + hp.elev)) * 0.5)
        r.tw = F32((float(self.tw + hp.tw)) * 0.5)
        r.wp = F32((float(self.wp + hp.wp)) * 0.5)
        r.cx = F32((float(self.cx + hp.cx)) * 0.5)
        r.n = F32((float(self.n + hp.n)) * 0.5)
        if self.elev < hp.elev:
            dz = F32(r.elev - self.elev)
            r.area = F32(float(self.area) + float(self.tw + r.tw) * 0.5 * float(dz))
        else:
            dz = F32(r.elev - hp.elev)
            r.area = F32(float(hp.area) + float(hp.tw + r.tw) * 0.5 * float(dz))
        return r


def _sample_elevations(z: np.ndarray, laminar: float, elev_tol: float) -> list:
    zs = sorted(z.tolist())
    zs = [F32(v) for v in zs]
    out = [zs[0]]
    lam = F32(laminar)
    if lam > 0:
        out.append(F32(zs[0] + lam))
    step = lam if float(lam) >= 0.001 else F32(0.001)
    for v in zs[1:]:
        if v > F32(out[-1] + step):
            out.append(v)
            out.append(F32(v + step))
    if out[-1] < zs[-1]:
        out.append(F32(zs[-1] + F32(elev_tol)))
    return out


def _profile_properties(s: np.ndarray, z: np.ndarray, elevs: list) -> list:
    """ComputeProfileHydraulicProperties: area, top width, wetted perimeter,
    centroid and ground intersections at each elevation (float32)."""
    E = np.array(elevs, dtype=F32)
    n = len(E)
    area = np.zeros(n, F32)
    tw = np.zeros(n, F32)
    wp = np.zeros(n, F32)
    cxs = np.zeros(n, F32)
    gri = np.zeros(n, np.int64)
    first = np.arange(n) == 0
    for j in range(1, len(s)):
        x, x2 = s[j - 1], s[j]
        y_a, y_b = z[j - 1], z[j]
        w = F32(x2 - x)
        mid = F32(float(x) + float(w) * 0.5)
        if y_a == y_b:
            sel = ~(E < y_a) & ((E > y_a) | first)
            a_ = (E - y_a) * w
            area[sel] += a_[sel]
            tw[sel] += w
            wp[sel] += w
            cxs[sel] += (a_ * mid)[sel]
            continue
        rising = y_a < y_b
        lo, hi = (y_a, y_b) if rising else (y_b, y_a)
        seg_len = F32(np.sqrt(float(x2 - x) ** 2 + float(y_b - y_a) ** 2))
        rise = F32(hi - lo)
        tri = F32(F32(w * F32(hi - lo)) * F32(0.5))
        above = E >= hi
        part = (E > lo) & ~above
        a_full = (E - hi) * w
        area[above] += (tri + a_full)[above]
        tw[above] += w
        wp[above] += seg_len
        if rising:
            cxs[above] += (a_full * mid + tri * F32(x + w / F32(3)))[above]
        else:
            cxs[above] += (a_full * mid + tri * F32(x2 - w / F32(3)))[above]
        with np.errstate(divide="ignore", invalid="ignore"):
            dz = E - lo
            frac = dz / rise
            ww = w * frac
            a_p = ww * dz * F32(0.5)
        area[part] += a_p[part]
        tw[part] += ww[part]
        wp[part] += (frac * seg_len)[part]
        if rising:
            cxs[part] += (a_p * (x + ww / F32(3)))[part]
        else:
            cxs[part] += (a_p * (x2 - ww / F32(3)))[part]
        gri[(E > lo) & (E <= hi)] += 1
    hps = [_HP(E[k], area[k], tw[k], wp[k], cxs[k], 0, int(gri[k])) for k in range(n)]
    last = F32(0)
    for hp in reversed(hps):
        if hp.area > 0:
            hp.cx = F32(hp.cx / hp.area)
            last = hp.cx
        else:
            hp.cx = last
    return hps


def _ensure_tolerance(elev_tol, k_tol, lo, hi, out):
    k_lo, k_hi = lo.conveyance(), hi.conveyance()
    if (k_lo < k_hi) and (F32(hi.elev - lo.elev) > F32(elev_tol)):
        base = F32((float(k_lo + k_hi)) * 0.5)
        mid = lo.midpoint(hi)
        if F32(abs(base - mid.conveyance()) / base) > F32(k_tol):
            _ensure_tolerance(elev_tol, k_tol, lo, mid, out)
    out.append(hi)


def _interp_conveyance(elev, below, above):
    dz = F32(F32(elev) - below.elev)
    af = F32(dz / F32(above.elev - below.elev))
    bf = F32(1.0 - float(af))
    return F32(below.conveyance() * bf + above.conveyance() * af)


def _filter_on_conveyance(k_tol, hps):
    keep = {0, len(hps) - 1}
    stack = [(0, len(hps) - 1)]
    while stack:
        first, last = stack.pop()
        if last - first <= 1:
            continue
        best, idx = F32(0), 0
        for i in range(first + 1, last):
            k = hps[i].conveyance()
            ck = _interp_conveyance(hps[i].elev, hps[first], hps[last])
            r = F32(abs(k - ck) / k)
            if r > best:
                best, idx = r, i
        if best > F32(k_tol):
            keep.add(idx)
            stack.append((idx, last))
            stack.append((first, idx))
    return [hps[i] for i in sorted(keep)]


def _bisect_conveyance(lo, hi, ref, elev_tol):
    mid = lo.midpoint(hi)
    it = 1
    while not (F32(hi.elev - lo.elev) < F32(elev_tol)):
        if mid.conveyance() < ref.projected_conveyance(mid.elev):
            lo = mid
        else:
            hi = mid
        mid = lo.midpoint(hi)
        it += 1
        if it > 8:
            break
    return mid


def face_table(s, z, mann_n, elev_tol, k_tol, laminar) -> list:
    """ProfileHydraulicPropertiesTable (one Manning's n): the face's
    (elevation, area, wetted perimeter, centroid) rows."""
    zmin, zmax = z.min(), z.max()
    width = F32(s[-1] - s[0])
    n = F32(mann_n)
    if F32(zmax - zmin) <= 0:
        rows = [_HP(zmin, 0, width, width, F32(float(width) * 0.5), n)]
        if F32(laminar) != 0:
            rows.append(_HP(F32(zmin + F32(laminar)), F32(width * F32(laminar)), width, width,
                            F32(float(width) * 0.5), n))
        return rows
    geom = _profile_properties(s, z, _sample_elevations(z, laminar, elev_tol))
    for hp in geom:
        hp.n = n
    if k_tol <= 0:
        return geom
    return _refine(geom, width, laminar, elev_tol, k_tol)


def face_table_n_profile(s, z, n_stations, n_values, mesh_n, elev_tol, k_tol, laminar) -> list:
    """ProfileHydraulicPropertiesTable with a Manning's n profile along the
    face (Multiple Face Mann n): the face is split where n changes and each
    row's n is the wetted-perimeter composite of the parts. A row with n 0
    (the face didn't split) takes the face's own n when written."""
    zmin, zmax = z.min(), z.max()
    width = F32(s[-1] - s[0])
    subs, n_lo, n_hi = _n_subsections(s, z, n_stations, n_values, mesh_n, width)
    if F32(zmax - zmin) <= 0:
        n = F32(0)
        if subs:
            total_len, total = F32(0), 0.0
            for ps, _, sn in subs:
                length = F32(ps[-1] - ps[0])
                total += float(length) * float(sn) ** 1.5
                total_len = F32(total_len + length)
            n = F32((total / float(total_len)) ** (2.0 / 3.0))
        rows = [_HP(zmin, 0, width, width, F32(float(width) * 0.5), n)]
        if F32(laminar) > 0:
            rows.append(_HP(F32(zmin + F32(laminar)), F32(width * F32(laminar)), width, width,
                            F32(width * F32(0.5)), n))
        return rows
    elevs = _sample_elevations(z, laminar, elev_tol)
    if subs:
        parts = [(_profile_properties(ps, pz, elevs), sn) for ps, pz, sn in subs]
        geom = [None] * len(elevs)
        cx, n = F32(0), F32(0)
        for r in range(len(elevs) - 1, -1, -1):
            hp = _HP(elevs[r])
            moment, k = 0.0, 0.0
            for g, sn in parts:
                e = g[r]
                hp.area = F32(hp.area + e.area)
                hp.tw = F32(hp.tw + e.tw)
                hp.wp = F32(hp.wp + e.wp)
                moment += float(F32(e.area * e.cx))
                k += float(e.wp) * float(sn) ** 1.5
            if float(hp.area) > 0:
                cx = F32(moment / float(hp.area))
                n = F32((k / float(hp.wp)) ** 0.6666666865348816)
                n = n_lo if n < n_lo else (n_hi if n > n_hi else n)
            hp.cx, hp.n = cx, n
            geom[r] = hp
    else:
        geom = _profile_properties(s, z, elevs)
    if k_tol <= 0 or np.isnan(k_tol):
        return geom
    return _refine(geom, width, laminar, elev_tol, k_tol, top=zmax)


def _n_subsections(s, z, n_stations, n_values, mesh_n, width):
    """The face profile's parts, one per Manning's n reach, as (stations,
    elevations, n); None unless there are at least two. Also the lowest and
    highest n, which bound the composite."""
    xs = [F32(x) for x in n_stations]
    ns = [F32(mesh_n) if (np.isnan(v) or v <= 0) else F32(v) for v in n_values]
    if len(xs) <= 1:
        return None, None, None
    last = F32(s[-1])
    short = F32(width / F32(1000))
    for j in range(len(xs) - 1, -1, -1):
        nxt = last if j >= len(xs) - 1 else xs[j + 1]
        if F32(nxt - xs[j]) < short:
            del xs[j], ns[j]
    if len(xs) <= 1:
        return None, None, None
    xs[0] = F32(0)
    subs = []
    for k in range(len(xs)):
        nxt = last if k >= len(xs) - 1 else xs[k + 1]
        ps, pz = _sub_section(s, z, xs[k], nxt)
        if len(ps) <= 1 or not (ps[-1] > ps[0]):
            continue
        subs.append((np.array(ps, F32), np.array(pz, F32), ns[k]))
    if len(subs) <= 1:
        return None, None, None
    return subs, min(ns), max(ns)


def _sub_section(s, z, left, right):
    """Profile.GetSubSection: the profile between two stations, with points
    interpolated where it is cut."""
    ps, pz = [], []

    def add(x, y):
        if not ps or not (ps[-1] == x and pz[-1] == y):
            ps.append(x)
            pz.append(y)

    for i in range(len(s) - 1):
        xl, yl, xr, yr = s[i], z[i], s[i + 1], z[i + 1]
        if xl > right or xr < left:
            continue
        if xl == xr:
            if xl == left:
                ps.clear()
                pz.clear()
                if yl > yr:
                    add(xl, yl)
                    add(xr, yr)
            elif xl == right:
                if yl < yr:
                    add(xl, yl)
                    add(xr, yr)
                break
            else:
                add(xl, yl)
                add(xr, yr)
        elif xl >= left and xr <= right:
            add(xl, yl)
            add(xr, yr)
        else:
            m = F32(F32(yr - yl) / F32(xr - xl))
            if xl < left:
                add(F32(left), F32(yl + F32(m * F32(left - xl))))
            else:
                add(xl, yl)
            if xr > right:
                add(F32(right), F32(yl + F32(m * F32(right - xl))))
            else:
                add(xr, yr)
    return ps, pz


def _refine(geom, width, laminar, elev_tol, k_tol, top=None) -> list:
    """The rows kept from the sampled geometry: interpolation tolerance on
    conveyance, conveyance that must rise, then thinning; extended upward
    when the last row falls short of the full width (or, with an n profile,
    of the profile's top)."""
    lst = [geom[0]]
    start = 1
    if F32(laminar) > 0:
        lst.append(geom[1])
        start = 2
    for j in range(start, len(geom)):
        _ensure_tolerance(elev_tol, k_tol, geom[j - 1], geom[j], lst)
    out = [lst[0]]
    for i in range(1, len(lst)):
        item = lst[i]
        prev = out[-1]
        if not (item.conveyance() > prev.projected_conveyance(item.elev)):
            continue
        if prev.elev != lst[i - 1].elev:
            out.append(_bisect_conveyance(lst[i - 1], lst[i], prev, elev_tol))
        out.append(item)
    out = _filter_on_conveyance(k_tol, out)
    if out[-1].tw < lst[-1].tw or (top is not None and out[-1].elev < top):
        hi = lst[-1].copy()
        while True:
            lo = hi.copy()
            hi.elev = F32(hi.elev + F32(1))
            hi.area = F32(hi.area + width)
            if hi.conveyance() > out[-1].projected_conveyance(hi.elev):
                break
        out.append(_bisect_conveyance(lo, hi, out[-1], elev_tol))
    return out


def centroid_at_area(rows, fraction=0.05):
    """ComputeCentroidAtArea: centroid station of the bottom fraction of area."""
    target = F32(rows[-1].area * F32(fraction))
    if target <= rows[0].area:
        return rows[0].cx
    areas = [r.area for r in rows]
    lo, hi = 0, len(areas) - 1
    while lo <= hi:   # Array.BinarySearch
        mid = (lo + hi) // 2
        if areas[mid] == target:
            return rows[mid].cx
        if areas[mid] < target:
            lo = mid + 1
        else:
            hi = mid - 1
    i = lo
    f = F32(F32(target - areas[i - 1]) / F32(areas[i] - areas[i - 1]))
    return F32(rows[i].cx * f + rows[i - 1].cx * F32(F32(1) - f))


# ----------------------------------------------------------------------------
# Manning's n along faces (Multiple Face Mann n)
# ----------------------------------------------------------------------------

N_NODATA = F32(-9999)


def _seg_length(ax, ay, bx, by):
    """PointM.DistanceTo."""
    return math.pow(math.pow(ax - bx, 2.0) + math.pow(ay - by, 2.0), 0.5)


@dataclass
class ClassTiles:
    """A classification raster's tile grid as RASMapper reads it
    (LandCoverLayer.TileRaster): tiles from the GeoTIFF's origin, one more
    row and column than the file has, each read as 256 x 256 cells."""
    min_x: float
    max_y: float
    tile: float
    rows: int
    cols: int
    classes: object     # classes(min_x, max_y, cell_size) -> 256 x 256 class indexes

    @classmethod
    def from_tif(cls, tif: Path, classes) -> "ClassTiles":
        ds = gdal.Open(str(tif))
        gt = ds.GetGeoTransform()
        bw, bh = ds.GetRasterBand(1).GetBlockSize()
        w, h = ds.RasterXSize, ds.RasterYSize
        ds = None
        pixel = max(w * gt[1] / w, h * -gt[5] / h)
        span = bw * pixel
        ntw, ntt = -(-w // bw), -(-h // bh)
        tile = max(((gt[0] + ntw * span) - gt[0]) / ntw, (gt[3] - (gt[3] - ntt * span)) / ntt)
        return cls(gt[0], gt[3], tile, ntt + 1, ntw + 1, classes)


def face_n_profiles(polylines: list, tiles: ClassTiles, n_of_class, replace_nodata: float,
                    decimal_places: int = 1) -> list:
    """RASGeometry.ExtractManningsNPolylines, from the land cover raster (no
    regions): for each face polyline, (stations, n) where n changes, float32.
    n_of_class[c] is class c's n, -9999 where it has none; -9999 in the result
    becomes replace_nodata."""
    lengths = [[_seg_length(*p[j], *p[j + 1]) for j in range(len(p) - 1)] for p in polylines]
    raster = _n_from_raster(polylines, lengths, tiles, n_of_class)
    out = []
    for i, p in enumerate(polylines):
        total = 0.0
        for length in lengths[i]:
            total += length
        values = raster.get(i) or [(0.0, total, -1, -1, N_NODATA)]
        values = _fill_n_gaps(values, total, 1.0 / math.pow(10.0, decimal_places))
        out.append(_n_breaks(values, F32(replace_nodata)))
    return out


def _n_from_raster(polylines, lengths, tiles, n_of_class) -> dict:
    """ExtractLCVariablePolylinesRegionsRasterExtraction: per face, runs of
    one class as (start, end, class, -1, n)."""
    t = tiles
    max_x = t.min_x + t.tile * t.cols
    min_y = t.max_y - t.tile * t.rows
    work = {}     # tile -> [(face, segment, start station)], in the order first reached
    for i, p in enumerate(polylines):
        start = 0.0
        for j in range(len(p) - 1):
            (ax, ay), (bx, by) = p[j], p[j + 1]
            clipped = _clip(ax, ay, bx, by, t.min_x, max_x, min_y, t.max_y)
            if clipped is not None:
                for rc in _walk_cells(*clipped, t.min_x, max_x, min_y, t.max_y,
                                      t.tile, t.rows, t.cols):
                    work.setdefault(rc, []).append((i, j, start))
            start += lengths[i][j]
    hits = {}
    for (tr, tc), segs in work.items():
        x0, x1 = t.min_x + t.tile * tc, t.min_x + t.tile * (tc + 1)
        y1, y0 = t.max_y - t.tile * tr, t.max_y - t.tile * (tr + 1)
        cs = max((x1 - x0) / 256, (y1 - y0) / 256)
        classes = t.classes(x0, y1, cs)
        for i, j, start in segs:
            (ax, ay), (bx, by) = polylines[i][j], polylines[i][j + 1]
            length = lengths[i][j]
            for r, c in _walk_cells(ax, ay, bx, by, x0, x1, y0, y1, cs, 256, 256):
                if not (0 <= r < 256 and 0 <= c < 256):
                    continue
                lc = int(classes[r, c])
                px0, px1 = x0 + cs * c, x0 + cs * (c + 1)
                py1, py0 = y1 - cs * r, y1 - cs * (r + 1)
                corners = ((px0, py0), (px0, py1), (px1, py1), (px1, py0))
                found = []
                for k in range(4):
                    hit, _, _, t1, _ = tin._seg_intersect((ax, ay), (bx, by), corners[k],
                                                          corners[(k + 1) % 4])
                    if hit:
                        found.append(start + length * t1)
                inside = lambda x, y: px0 <= x <= px1 and py0 <= y <= py1
                if found:
                    a = found[0]
                    b = found[1] if len(found) > 1 else (start if inside(ax, ay) else start + length)
                    entry = (min(a, b), max(a, b))
                elif inside(ax, ay) and inside(bx, by):
                    entry = (start, start + length)
                else:
                    continue
                hits.setdefault(i, []).append((*entry, lc, -1, F32(n_of_class[lc])))
    runs = {}
    for i in range(len(polylines)):
        if i not in hits:
            continue
        lst = sorted(hits[i], key=lambda v: v[0])
        merged = [[0.0, lst[0][1], lst[0][2], -1, lst[0][4]]]
        for k, item in enumerate(lst):
            value = list(merged[-1])
            idx = len(merged) - 1
            if item[2] != value[2]:
                merged.append(list(item))
                value[1] = item[0]
                merged[idx] = value
                value = list(item)
                idx = len(merged) - 1
            if k == len(lst) - 1:
                value[1] = item[1]
                merged[idx] = value
        runs[i] = [tuple(v) for v in merged]
    return runs


def _fill_n_gaps(values, total, tol) -> list:
    """ExtractLCVariablePolylineFillNoData: no-data runs where the face has none."""
    out = list(values)
    end = 0.0
    lst = sorted(values, key=lambda v: v[0])
    for v in lst:
        if v[0] > end and v[0] - end > tol:
            out.append((end, v[0], -1, -1, N_NODATA))
        end = v[1]
    if not lst or end < total:
        out.append((end, total, -1, -1, N_NODATA))
    return out


def _n_breaks(values, replace_nodata) -> tuple:
    """ExtractLCVariablePolyineSort: the stations where n changes."""
    lst = sorted(values, key=lambda v: v[0])
    first = lst[0]
    xs, ns = [first[0]], [first[4]]
    prev_n, end, end_n = first[4], first[1], first[4]
    override, override_end = first[3] != -1, first[1] if first[3] != -1 else 0.0
    for v in lst:
        if v == first:
            continue
        if override and v[0] > override_end:
            xs.append(override_end)
            ns.append(end_n)
        if v[3] != -1:
            xs.append(v[0])
            ns.append(v[4])
            override_end, override = v[1], True
        elif v[4] != prev_n:
            xs.append(v[0])
            ns.append(v[4])
            override = False
        prev_n = v[4]
        if v[1] > end:
            end, end_n = v[1], v[4]
    s = np.array(xs, dtype=F32)
    n = np.array([replace_nodata if v == N_NODATA else v for v in ns], dtype=F32)
    return s, n


# ----------------------------------------------------------------------------
# Cells
# ----------------------------------------------------------------------------

def polygon_area(poly: np.ndarray, signed: bool = False) -> float:
    """Polygon.Area: trapezoids measured from the lowest y, in ring order."""
    ring = np.vstack([poly, poly[:1]]) if not np.array_equal(poly[0], poly[-1]) else poly
    min_y = ring[:, 1].min()
    total = 0.0
    for (ax, ay), (bx, by) in zip(ring[:-1], ring[1:]):
        total += ((ay - min_y) + (by - min_y)) * (bx - ax)
    return total * 0.5 if signed else abs(total * 0.5)


def _cell_raster(max_x, min_x, max_y, min_y, max_cells=CELL_RASTER_SIZE):
    """RasterM.ComputeRasterM: (rows, cols, min_x, max_y, cell size)."""
    w, h = max_x - min_x, max_y - min_y
    ratio = w / h
    k = np.sqrt(max_cells / ratio)
    cols = int(np.floor(k * ratio))
    rows = int(np.floor(k))
    if rows == 0:
        rows, cols = 1, max_cells
    elif cols == 0:
        cols, rows = 1, max_cells
    cs = max(h / rows, w / cols) * 1.01
    cx, cy = min_x + 0.5 * w, min_y + 0.5 * h
    hx, hy = 0.5 * cols * cs, 0.5 * rows * cs
    max_x, min_x, max_y, min_y = cx + hx, cx - hx, cy + hy, cy - hy
    cs = max((max_x - min_x) / cols, (max_y - min_y) / rows)
    return rows, cols, min_x, max_y, max_x, min_y, cs


def _rasterize(poly: np.ndarray, rows, cols, min_x, max_y, max_x, cs) -> np.ndarray:
    """RasterizePolygon.ComputeCells: raster cells whose centre is inside."""
    mask = np.zeros((rows, cols), bool)
    ring = np.vstack([poly, poly[:1]])
    yc = max_y - cs * (np.arange(rows) + 0.5)
    xs = [[] for _ in range(rows)]
    eps = cs * 1e-5
    centre = lambda y: max_y - cs * (np.floor((max_y - y) / cs) + 0.5)
    for (x1, y1), (x2, y2) in zip(ring[:-1], ring[1:]):
        if y1 == y2:
            continue
        lo_y, hi_y = min(y1, y2), max(y1, y2)
        lo_x, hi_x = (x1, x2) if y1 < y2 else (x2, x1)
        slope = (hi_x - lo_x) / (hi_y - lo_y)
        # An end exactly on a row centre is moved down a little, so each edge
        # counts at its lower end and not at its upper end.
        if hi_y == centre(hi_y):
            hi_y -= eps
        if lo_y == centre(lo_y):
            lo_y -= eps
        inside = np.where((yc > lo_y) & (yc < hi_y))[0]
        for r in inside:
            x = lo_x + (yc[r] - lo_y) * slope
            if x < max_x:
                xs[r].append(x)
    xcen = min_x + cs * (np.arange(cols) + 0.5)
    for r in range(rows):
        row = sorted(xs[r])
        if len(row) % 2:
            row.append(max_x)
        for x0, x1 in zip(row[0::2], row[1::2]):
            mask[r] |= (xcen >= x0) & (xcen <= x1)
    return mask


def cell_table(terrain: Terrain, poly: np.ndarray, area: float, lowest_face: float,
               vol_tol: float, min_area_frac: float) -> list:
    """The volume-elevation table of one cell: [(elevation, volume), ...]."""
    min_x, min_y = poly.min(axis=0)
    max_x, max_y = poly.max(axis=0)
    rows, cols, rmin_x, rmax_y, rmax_x, rmin_y, cs = _cell_raster(max_x, min_x, max_y, min_y)
    xc = rmin_x + cs * (np.arange(cols) + 0.5)
    yc = rmax_y - cs * (np.arange(rows) + 0.5)
    X, Y = np.meshgrid(xc, yc)
    # TerrainLayer.Resample: each raster fills what the ones before it left
    # empty, at its own pyramid level; stitches from the render blocks under
    # the raster, two blocks either side, are bridged by a TIN painted over.
    vals = np.full(rows * cols, np.nan, dtype=F32)
    segs = []
    for raster in _as_layer(terrain).rasters:
        level = raster.level_for(cs)
        v = raster.sample(X.ravel(), Y.ravel(), level)
        vals = np.where(np.isnan(vals), v, vals).astype(F32)
        lv = raster.levels[level]
        bc = np.floor((xc[[0, -1]] - lv.x0) / lv.cell).astype(int)
        br = np.floor((lv.y0 - yc[[0, -1]]) / lv.cell).astype(int)
        found = _stitch_segments(lv, br.min() - 2, br.max() + 2, bc.min() - 2, bc.max() + 2)
        # ResampleInternal orders each raster's stitches by SegmentM.GetHashCode
        segs.extend(sorted(found, key=lambda sg: tin.segment_hash(sg[0])))
    if len(segs) > 1:
        vals = np.where(np.isfinite(vals), vals, F32(NODATA)).astype(F32)
        tin.paint(tin.StitchTin(segs).compute(), vals, rows, cols, rmin_x, rmax_y, cs, F32)
        vals = np.where(vals == F32(NODATA), np.nan, vals).astype(F32)
    inside = _rasterize(poly, rows, cols, rmin_x, rmax_y, rmax_x, cs).ravel()
    sel = inside & np.isfinite(vals)
    if not sel.any():
        return []
    v = vals[sel]
    vmin, vmax = v.min(), v.max()
    if vmax > vmin:
        # Histogram of wet area by elevation, 256 bins, values added in runs.
        inc = (float(vmax) - float(vmin)) / HISTOGRAM_BINS
        off = float(vmin) + inc * 0.5
        vd = v.astype(np.float64)
        idx = np.rint((vd - off) * (1.0 / inc)).astype(np.int64)
        idx[vd <= float(vmin)] = 0
        idx[vd >= float(vmax)] = HISTOGRAM_BINS - 1
        bins = np.zeros(HISTOGRAM_BINS)
        pix = cs * cs
        # run-length grouping of equal consecutive values, as RASMapper adds them
        starts = np.r_[0, np.where(v[1:] != v[:-1])[0] + 1]
        counts = np.diff(np.r_[starts, len(v)])
        np.add.at(bins, idx[starts], pix * counts)
        cum = np.cumsum(bins).astype(F32)               # area up to each bin
        elev = (float(vmin) + inc * np.arange(HISTOGRAM_BINS)).astype(F32)  # bin minimum
        ratio = F32(area / float(cum[-1]))
        cum = (cum * ratio).astype(F32)
        mid = ((cum[:-1] + cum[1:]).astype(np.float64) * 0.5).astype(F32)
        dz = (elev[1:] - elev[:-1]).astype(F32)
        vol = np.add.accumulate((dz * mid).astype(F32), dtype=F32)
        px = np.r_[F32(0), vol].astype(F32)
        py = np.r_[F32(vmin), elev[1:]].astype(F32)
        if vol_tol > 0:
            idx = _dp_keep(px, py, vol_tol)
            px, py = px[idx], py[idx]
        pts = [[px[i], py[i]] for i in range(len(px))]
        if min_area_frac > 0:
            pts = _apply_min_area_fraction(pts, area, min_area_frac, vol_tol)
    else:
        pts = [[F32(0), F32(vmin)]]
    if F32(lowest_face) < pts[0][1]:
        pts[0][1] = F32(lowest_face)
    return [(float(p[1]), float(p[0])) for p in pts]


def _apply_min_area_fraction(pts, area, frac, vol_tol):
    """Raise the curve's start so the bottom is where frac of the cell is wet."""
    target = area * float(F32(frac))
    for i in range(1, len(pts)):
        if not (float(pts[i][0] / F32(pts[i][1] - pts[0][1])) > target):
            continue
        if i == 1:
            break
        ax, ay = float(pts[i - 1][0]), float(pts[i - 1][1])
        bx, by = float(pts[i][0]), float(pts[i][1])
        cx, cy = 0.0, float(pts[0][1])
        dx, dy = float(F32(pts[i][1] - pts[0][1])) * target, float(pts[i][1])
        hit, t1, _ = _segment_intersections(ax, ay, bx, by, np.array([cx]), np.array([cy]),
                                            np.array([dx]), np.array([dy]))
        if hit[0]:
            x = ax + t1[0] * (bx - ax)
            y = ay + t1[0] * (by - ay)
            if i > 2:
                del pts[1:i - 1]
            if float(pts[2][1]) - y > float(F32(vol_tol)):
                pts[1] = [F32(x), F32(y)]
            else:
                del pts[1]
        break
    return pts
