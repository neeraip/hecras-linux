"""HEC-RAS 2D geometry preprocessor for Linux.

Prepares a project for HEC's Linux engines without the Windows GUI:

  A  A GUI-computed project: copies the geometry, strips the old results.
  B  A RASMapper mesh without property tables: computes them exactly as
     RASMapper does (property_tables.py), with land cover and infiltration.

Input:  .g01.hdf + .p01 + terrain + Land Classification/
Output: p01.tmp.hdf ready for RasGeomPreprocess -> RasUnsteady

Usage:
    hecras-preprocess /path/to/project_dir [--project-name Muncie] [--plan 01]
    python ras_preprocess.py ...   (the same, from a clone)
"""

import argparse
import logging
import shutil
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import h5py
import numpy as np
from osgeo import gdal

from . import property_tables

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Suppress GDAL warnings
gdal.UseExceptions()

# ============================================================================
# DATA CLASSES
# ============================================================================

@dataclass
class MeshData:
    """A 2D flow area's mesh, as read from a RASMapper geometry HDF."""
    cell_centers: np.ndarray = None      # (n_cells, 2) all cells
    cell_points_input: np.ndarray = None # (n_interior, 2) user input points
    face_points: np.ndarray = None       # (n_fp, 2) vertices
    perimeter: np.ndarray = None         # (n_perim, 2)

    # Connectivity
    faces_cell_idx: np.ndarray = None    # (n_faces, 2)
    faces_fp_idx: np.ndarray = None      # (n_faces, 2)
    cells_fp_idx: np.ndarray = None      # (n_cells, max_fp) padded with -1
    fp_is_perimeter: np.ndarray = None   # (n_fp,)

    # Derived
    faces_normal_length: np.ndarray = None  # (n_faces, 3) [nx, ny, length]
    cell_is_boundary: np.ndarray = None     # (n_cells,) bool
    # Perimeter vertices lying between a perimeter face's two face points
    # (HEC-RAS "Faces Perimeter Info/Values"); one (k, 2) array per face.
    faces_perimeter_pts: list = None

    @property
    def n_cells(self): return len(self.cell_centers)
    @property
    def n_faces(self): return len(self.faces_cell_idx)
    @property
    def n_fp(self): return len(self.face_points)
    @property
    def n_interior(self): return int((~self.cell_is_boundary).sum())

# ============================================================================
# STEP 3: PROPERTY TABLES (property_tables.py)
# ============================================================================

def _face_polyline(mesh: MeshData, fi: int) -> np.ndarray:
    """A face's points: its two face points with any perimeter vertices between."""
    mid = mesh.faces_perimeter_pts[fi] if mesh.faces_perimeter_pts else []
    return np.vstack([mesh.face_points[mesh.faces_fp_idx[fi, 0]],
                      np.asarray(mid, dtype=np.float64).reshape(-1, 2),
                      mesh.face_points[mesh.faces_fp_idx[fi, 1]]])


def _cell_polygon(mesh: MeshData, ci: int, face_of_pair: dict) -> np.ndarray:
    """A cell's outline, including the perimeter vertices of its bent faces."""
    fps = mesh.cells_fp_idx[ci]
    fps = fps[fps >= 0]
    ring = []
    for a, b in zip(fps, np.roll(fps, -1)):
        ring.append(mesh.face_points[a])
        fi = face_of_pair.get((int(a), int(b)))
        if fi is not None and mesh.faces_perimeter_pts and len(mesh.faces_perimeter_pts[fi]):
            mid = np.asarray(mesh.faces_perimeter_pts[fi], dtype=np.float64).reshape(-1, 2)
            if mesh.faces_fp_idx[fi, 0] != a:
                mid = mid[::-1]
            ring.extend(mid)
    return np.array(ring)


def table_settings(g01_hdf_path: Path = None, area_name: str = "") -> dict:
    """The 2D area's property-table settings from the Attributes table of its
    g01.hdf, which RASMapper writes with the mesh; defaults otherwise."""
    settings = dict(property_tables.DEFAULTS, composite_lc=False)
    if g01_hdf_path is None:
        return settings
    fields = {"cell_volume_tolerance": "Cell Vol Tol",
              "cell_min_area_fraction": "Cell Min Area Fraction",
              "face_profile_tolerance": "Face Profile Tol",
              "face_area_elevation_tolerance": "Face Area Tol",
              "face_conveyance_ratio": "Face Conv Ratio",
              "laminar_depth": "Laminar Depth"}
    with h5py.File(str(g01_hdf_path), "r") as hf:
        rows = hf["Geometry/2D Flow Areas/Attributes"][()]
    for row in rows:
        if row["Name"].decode("utf-8", errors="replace").strip("\x00").strip() != area_name:
            continue
        for key, field in fields.items():
            if field in rows.dtype.names and np.isfinite(row[field]):
                settings[key] = float(row[field])
        settings["composite_lc"] = bool("Composite LC" in rows.dtype.names and row["Composite LC"])
        for flag in ("Multiple Face Mann n",):
            if flag in rows.dtype.names and row[flag]:
                raise SystemExit(f"ras_preprocess.py: '{flag}' is on for '{area_name}'; "
                                 "property tables for it aren't supported. "
                                 "Turn it off in RASMapper, or use Workflow A.")
    return settings


def compute_face_properties(mesh: MeshData, terrain: "property_tables.Terrain",
                            mannings_per_face: np.ndarray, settings: dict) -> dict:
    """Per-face minimum elevation, low-elevation centroid and (elevation, area,
    wetted perimeter, Manning's n) table, as RASMapper computes them."""
    logger.info("Computing face properties from terrain...")
    n_faces = mesh.n_faces
    min_elevs = np.full(n_faces, property_tables.NAN, dtype=np.float32)
    low_elev_centroid = np.zeros(n_faces, dtype=np.float32)
    face_ae_info = np.zeros((n_faces, 2), dtype=np.int32)
    face_ae_values = []

    for fi in range(n_faces):
        stations, elevs = property_tables.face_profile(
            terrain, _face_polyline(mesh, fi), settings["face_profile_tolerance"])
        if len(stations) == 0:
            face_ae_info[fi] = [len(face_ae_values), 0]
            continue
        n = np.float32(mannings_per_face[fi])
        rows = property_tables.face_table(stations, elevs, n,
                                          settings["face_area_elevation_tolerance"],
                                          settings["face_conveyance_ratio"],
                                          settings["laminar_depth"])
        min_elevs[fi] = rows[0].elev
        low_elev_centroid[fi] = property_tables.centroid_at_area(rows)
        face_ae_info[fi] = [len(face_ae_values), len(rows)]
        face_ae_values.extend([r.elev, r.area, r.wp, n] for r in rows)

    missing = int((face_ae_info[:, 1] == 0).sum())
    if missing:
        logger.warning("  %d face(s) have no terrain under them", missing)
    logger.info("  Computed %d faces, %d area-elevation points total",
                n_faces, len(face_ae_values))
    return {
        "min_elevations": min_elevs,
        "low_elev_centroid": low_elev_centroid,
        "face_ae_info": face_ae_info,
        "face_ae_values": np.array(face_ae_values, dtype=np.float32).reshape(-1, 4)
                          if face_ae_values else np.empty((0, 4), dtype=np.float32),
    }


# ============================================================================
# STEP 4: CELL VOLUME-ELEVATION TABLES
# ============================================================================

def compute_cell_properties(mesh: MeshData, terrain: "property_tables.Terrain",
                            face_min_elevs: np.ndarray, settings: dict) -> dict:
    """Per-cell minimum elevation, surface area and volume-elevation table, as
    RASMapper computes them. Needs the faces' minimum elevations: a cell's
    minimum is never above its lowest face."""
    logger.info("Computing cell properties from terrain...")
    n_cells = mesh.n_cells
    min_elevs = np.full(n_cells, property_tables.NAN, dtype=np.float32)
    surface_areas = np.zeros(n_cells, dtype=np.float32)
    vol_elev_info = np.zeros((n_cells, 2), dtype=np.int32)
    vol_elev_values = []

    face_of_pair = {}
    lowest_face = np.full(n_cells, np.inf, dtype=np.float32)
    for fi, (a, b) in enumerate(mesh.faces_fp_idx):
        face_of_pair[(int(a), int(b))] = face_of_pair[(int(b), int(a))] = fi
        for ci in mesh.faces_cell_idx[fi]:
            if 0 <= ci < n_cells and np.isfinite(face_min_elevs[fi]):
                lowest_face[ci] = min(lowest_face[ci], face_min_elevs[fi])

    for ci in range(n_cells):
        vol_elev_info[ci] = [len(vol_elev_values), 0]
        if mesh.cell_is_boundary[ci]:
            # RASMapper's outline of a ghost cell runs out along its one face
            # and back, so its area is zero up to rounding.
            fi = int(np.where((mesh.faces_cell_idx == ci).any(axis=1))[0][0])
            line = _face_polyline(mesh, fi)
            ring = np.vstack([line, line[1:-1][::-1]])
            surface_areas[ci] = property_tables.polygon_area(ring, signed=True)
            continue
        polygon = _cell_polygon(mesh, ci, face_of_pair)
        area = property_tables.polygon_area(polygon)
        surface_areas[ci] = area
        table = property_tables.cell_table(terrain, polygon, area, lowest_face[ci],
                                           settings["cell_volume_tolerance"],
                                           settings["cell_min_area_fraction"])
        if not table:
            continue
        min_elevs[ci] = table[0][0]
        vol_elev_info[ci] = [len(vol_elev_values), len(table)]
        vol_elev_values.extend(table)

    logger.info("  Computed %d cells, %d volume-elevation points total",
                n_cells, len(vol_elev_values))
    return {
        "min_elevations": min_elevs,
        "surface_areas": surface_areas,
        "vol_elev_info": vol_elev_info,
        "vol_elev_values": np.array(vol_elev_values, dtype=np.float32).reshape(-1, 2)
                           if vol_elev_values else np.empty((0, 2), dtype=np.float32),
    }



# ============================================================================
# STEP 5: LAND COVER SAMPLING
# ============================================================================

def _default_land_cover(mesh: MeshData, default_mannings: float) -> dict:
    """Land cover data when the project has none: one Manning's n everywhere."""
    return {"mannings_n": np.full(mesh.n_cells, default_mannings, dtype=np.float32),
            "face_mannings_n": np.full(mesh.n_faces, default_mannings, dtype=np.float32),
            "classifications": np.zeros(mesh.n_cells, dtype=np.int32),
            "face_classifications": np.zeros(mesh.n_faces, dtype=np.int32),
            "calibration_table": [("Default", default_mannings)],
            "id_to_name": {}, "name_to_imperv": {}, "id_to_mann": {}}


def sample_land_cover(mesh: MeshData, lc_hdf_path: str, lc_tif_path: str,
                      default_mannings: float = 0.06, composite: bool = False) -> dict:
    """Sample Manning's n and land cover properties per cell.

    LandCover.hdf structure (File Type = 'HEC Land Cover'):
      /Raster Map  - compound dataset: [('ID', i4), ('Name', S30)]
      /Variables   - compound dataset: [('Name', S30), ('ManningsN', f4), ('Percent Impervious', f4)]

    The raster pixel values correspond to IDs in /Raster Map.
    Manning's n is looked up via: pixel_id -> Name (via Raster Map) -> ManningsN (via Variables).
    """
    logger.info("Sampling land cover from: %s", lc_tif_path)

    ids = None
    names_rmap = None
    mann_vals = None
    imperv_vals = None
    var_names = None
    cal_table = []

    try:
        with h5py.File(lc_hdf_path, 'r') as hf:
            # /Raster Map: compound dataset at root level
            if 'Raster Map' in hf:
                rmap = hf['Raster Map'][()]          # structured numpy array
                ids = rmap['ID']                     # int32 array: raster pixel values
                names_rmap = rmap['Name']            # S30 array: class names

            # /Variables: compound dataset at root level
            if 'Variables' in hf:
                vdata = hf['Variables'][()]          # structured numpy array
                var_names = vdata['Name']            # S30: class names (same order as Variables)
                mann_vals = vdata['ManningsN']       # float32: Manning's n per class
                imperv_vals = vdata['Percent Impervious']  # float32: impervious fraction

        if ids is None or names_rmap is None or mann_vals is None:
            raise ValueError("Missing expected datasets in LandCover HDF")

        # Build name -> ManningsN lookup from /Variables
        name_to_mann = {}
        name_to_imperv = {}
        for i in range(len(var_names)):
            n_str = var_names[i].decode('utf-8', errors='replace').strip() if isinstance(var_names[i], bytes) else str(var_names[i]).strip()
            name_to_mann[n_str] = float(mann_vals[i])
            name_to_imperv[n_str] = float(imperv_vals[i])

        # Build ID -> Manning's n + calibration table from /Raster Map
        id_to_mann = {}
        id_to_name = {}
        for i in range(len(ids)):
            rmap_name = names_rmap[i].decode('utf-8', errors='replace').strip() if isinstance(names_rmap[i], bytes) else str(names_rmap[i]).strip()
            lc_id = int(ids[i])
            mann_n = name_to_mann.get(rmap_name, default_mannings)
            id_to_mann[lc_id] = mann_n
            id_to_name[lc_id] = rmap_name
            cal_table.append((rmap_name, mann_n))

        logger.info("  Loaded %d land cover classes from HDF", len(ids))

    except Exception as e:
        logger.warning("Could not read land cover HDF: %s. Using default Manning's n.", e)
        id_to_mann = {}
        id_to_name = {}
        cal_table = [("Default", default_mannings)]

    # Read land cover raster
    try:
        ds = gdal.Open(lc_tif_path)
        lc_raster = ds.GetRasterBand(1).ReadAsArray()
        lc_gt = ds.GetGeoTransform()
        ds = None
    except Exception as e:
        logger.warning("Could not read land cover TIF: %s. Using defaults.", e)
        return {**_default_land_cover(mesh, default_mannings), "calibration_table": cal_table}

    # Sample per cell center (interior cells only; boundary cells get default)
    n = mesh.n_cells
    mannings = np.full(n, default_mannings, dtype=np.float32)
    classifications = np.zeros(n, dtype=np.int32)

    raster_rows, raster_cols = lc_raster.shape

    # Every cell takes the land cover at its centre and every face at its
    # midpoint, as RASMapper does. Off the raster, or on a value the raster map
    # doesn't list, it is RASMapper's NoData class: the raster map's first row.
    nodata_id = int(ids[0]) if ids is not None and len(ids) else 0
    nodata_mann = id_to_mann.get(nodata_id, default_mannings)

    def lookup(points):
        cols = np.floor((points[:, 0] - lc_gt[0]) / lc_gt[1]).astype(int)
        rows = np.floor((points[:, 1] - lc_gt[3]) / lc_gt[5]).astype(int)
        inside = (rows >= 0) & (rows < raster_rows) & (cols >= 0) & (cols < raster_cols)
        lc_ids = np.full(len(points), nodata_id, dtype=np.int64)
        lc_ids[inside] = lc_raster[rows[inside], cols[inside]]
        lc_ids[~np.isin(lc_ids, list(id_to_mann))] = nodata_id
        return lc_ids, np.array([id_to_mann.get(int(i), nodata_mann) for i in lc_ids],
                                dtype=np.float32)

    cell_ids, mannings = lookup(mesh.cell_centers)
    classifications[~mesh.cell_is_boundary] = cell_ids[~mesh.cell_is_boundary]
    mid = (mesh.face_points[mesh.faces_fp_idx[:, 0]] + mesh.face_points[mesh.faces_fp_idx[:, 1]]) / 2
    face_ids, face_mannings = lookup(mid)

    percent_impervious = None
    if composite and ids is not None:
        # Interior cells average their whole raster instead (faces keep midpoints).
        layer = _ClassRaster(Path(lc_hdf_path))
        if Path(lc_tif_path) != layer.tif:
            layer.tif = Path(lc_tif_path)
            layer.raster, layer.gt = lc_raster, lc_gt
        names = [id_to_name.get(int(i), "") for i in ids]
        table = {"ManningsN": np.array([name_to_mann.get(n, np.nan) for n in names], dtype=np.float32),
                 "Percent Impervious": np.array([name_to_imperv.get(n, np.nan) for n in names], dtype=np.float32)}
        nd = names.index("NoData") if "NoData" in names else 0
        cls, avg = _composite_cells(mesh, {"lc": (layer.classes_at, table, nd, ids)})["lc"]
        mannings[:mesh.n_interior] = avg["ManningsN"]
        classifications[:mesh.n_interior] = cls
        percent_impervious = avg["Percent Impervious"]

    n_assigned = int(np.sum(mannings != default_mannings))
    logger.info("  Manning's n range: %.4f – %.4f  (%d/%d cells from land cover)",
                float(mannings.min()), float(mannings.max()), n_assigned, n)
    return {
        "mannings_n": mannings,
        "face_mannings_n": face_mannings,
        "classifications": classifications,
        "face_classifications": face_ids.astype(np.int32),
        "percent_impervious": percent_impervious,
        "calibration_table": cal_table,
        "id_to_name": id_to_name,
        "name_to_imperv": name_to_imperv if 'name_to_imperv' in dir() else {},
        "id_to_mann": id_to_mann,
    }


class _ClassRaster:
    """A RASMapper classification layer (LandCover.hdf, Soils.hdf, ...): its
    Raster Map IDs and raster, read once."""

    def __init__(self, layer_hdf: Path):
        with h5py.File(str(layer_hdf), "r") as hf:
            rmap = hf["Raster Map"][()]
        self.ids = rmap["ID"].astype(np.int64)
        self.names = [n.decode("utf-8", "replace").strip() for n in rmap["Name"]]
        self.tif = layer_hdf.with_suffix(".tif")
        ds = gdal.Open(str(self.tif))
        self.raster = ds.GetRasterBand(1).ReadAsArray()
        self.gt = ds.GetGeoTransform()
        ds = None
        self._levels = None
        # pixel value -> class index (position in the Raster Map); 0 (NoData) otherwise
        self._lookup = {int(v): i for i, v in enumerate(self.ids)}

    def _level_for(self, cell_size: float):
        """RasterPyramid.GetBestPyramidLevel: the finest level, unless the target
        pixel is at least the next level's size; returns (raster, geotransform)."""
        if self._levels is None:
            ds = gdal.Open(str(self.tif))
            band = ds.GetRasterBand(1)
            width = ds.RasterXSize * self.gt[1]
            self._levels = [(self.raster, self.gt)]
            for i in range(band.GetOverviewCount()):
                ov = band.GetOverview(i)
                cell = width / ov.XSize
                self._levels.append((ov.ReadAsArray(), (self.gt[0], cell, 0.0, self.gt[3], 0.0, -cell)))
            ds = None
        for i, (_, gt) in enumerate(self._levels):
            if gt[1] >= cell_size:
                return self._levels[max(i - 1, 0)]
        return self._levels[-1]

    def classes_at(self, points: np.ndarray, cell_size: float = None) -> np.ndarray:
        """Class index at each point: the pixel's ID's position in the Raster
        Map, 0 (the NoData class) off the raster or on an unlisted value. With
        cell_size, the pixel comes from the pyramid level RASMapper resamples a
        raster of that pixel size from."""
        r, gt = (self.raster, self.gt) if cell_size is None else self._level_for(cell_size)
        cols = np.floor((points[:, 0] - gt[0]) / gt[1]).astype(np.int64)
        rows = np.floor((points[:, 1] - gt[3]) / gt[5]).astype(np.int64)
        inside = (rows >= 0) & (rows < r.shape[0]) & (cols >= 0) & (cols < r.shape[1])
        classes = np.zeros(len(points), dtype=np.int64)
        vals = r[rows[inside], cols[inside]].astype(np.int64)
        uniq, inv = np.unique(vals, return_inverse=True)
        classes[inside] = np.array([self._lookup.get(int(v), 0) for v in uniq], dtype=np.int64)[inv]
        return classes


def _classify_points(layer_hdf: Path, points: np.ndarray) -> tuple:
    """A RASMapper classification layer's class at each point (see _ClassRaster).
    Returns (class indexes, number of classes)."""
    layer = _ClassRaster(layer_hdf)
    return layer.classes_at(points), len(layer.ids)


def _composite_cells(mesh: MeshData, samplers: dict) -> dict:
    """RASMapper's cell-averaged land cover (Composite LC, the CreatePropertyTables
    loop over ComputeAssociationCellAverageValue): each layer is sampled at the
    pixel centres of the cell's raster (the same as its volume table's), its
    variables averaged over the cell's pixels and its class the most common.

    RASMapper never clears the pixel mask between cells: each thread's mask
    keeps every pixel of the cells it handled before, and Parallel.SmartFor
    hands out cells in chunks of 1,000. That is reproduced (a fresh mask per
    chunk, assuming one thread per chunk), so values match RASMapper's; on a
    mesh of more than 1,000 cells RASMapper's own result depends on how its
    threads pick up the chunks.

    samplers: name -> (classes(points) -> class indexes, {variable: value per
    class index}, NoData class index, class IDs). Returns name -> (class ID per
    interior cell, {variable: float32 average per interior cell})."""
    out = {name: (np.zeros(mesh.n_interior, dtype=np.int32),
                  {v: np.zeros(mesh.n_interior, dtype=np.float32) for v in table})
           for name, (_, table, _, _) in samplers.items()}
    face_of_pair = {}
    for fi, (a, b) in enumerate(mesh.faces_fp_idx):
        face_of_pair[(int(a), int(b))] = face_of_pair[(int(b), int(a))] = fi
    mask = np.zeros(property_tables.CELL_RASTER_SIZE, dtype=bool)
    for ci in range(mesh.n_interior):
        if ci % 1000 == 0:
            mask[:] = False
        polygon = _cell_polygon(mesh, ci, face_of_pair)
        x0, y0 = polygon.min(axis=0)
        x1, y1 = polygon.max(axis=0)
        rows, cols, min_x, max_y, max_x, _, cs = property_tables._cell_raster(x1, x0, y1, y0)
        n = rows * cols
        mask[:n] |= property_tables._rasterize(polygon, rows, cols, min_x, max_y, max_x, cs).ravel()
        cell_mask = mask[:n]
        xc = min_x + cs * (np.arange(cols) + 0.5)
        yc = max_y - cs * (np.arange(rows) + 0.5)
        X, Y = np.meshgrid(xc, yc)
        points = np.column_stack([X.ravel(), Y.ravel()])
        for name, (classes_at, table, nd, class_ids) in samplers.items():
            classes = classes_at(points, cs)
            cells_cls, cells_var = out[name]
            for var, per_class in table.items():
                vals = per_class[classes][cell_mask]
                keep = ~np.isnan(vals) & (vals != np.float32(-9999)) & (vals != np.finfo(np.float32).max)
                total = np.add.reduce(vals[keep], dtype=np.float32) if keep.sum() < 2 else \
                    np.add.accumulate(vals[keep], dtype=np.float32)[-1]
                with np.errstate(invalid="ignore", divide="ignore"):
                    cells_var[var][ci] = np.float32(total) / np.float32(keep.sum())
            c = classes[cell_mask]
            counts = np.bincount(c[c != nd], minlength=len(class_ids))
            cells_cls[ci] = class_ids[int(np.argmax(counts))]
    return out


def sample_infiltration(mesh: MeshData, inf_hdf_path: str,
                        lc_tif_path: str = None, composite: bool = False) -> dict:
    """Infiltration parameters per interior cell and class per face, as
    RASMapper assigns them: the infiltration class at a cell's centre (a face's
    midpoint) selects a row of the layer's Variables table.

    An infiltration layer is either its own raster or "virtual": a combination
    of other layers (land cover and soils), whose class is the dependencies'
    class indexes combined in order, each scaled by the class counts after it.
    """
    logger.info("Sampling infiltration properties...")
    result = {"available": False, "model_type": "", "variables": {},
              "cell_classifications": np.zeros(mesh.n_interior, dtype=np.int32),
              "face_classifications": np.zeros(mesh.n_faces, dtype=np.int32),
              "properties": None}
    if not inf_hdf_path:
        return result
    inf_hdf = Path(inf_hdf_path)
    with h5py.File(str(inf_hdf), "r") as hf:
        result["model_type"] = hf.attrs.get("LC Type", b"").decode("utf-8", errors="replace").strip()
        deps = [d.decode("utf-8", errors="replace") for d in hf.attrs.get("Virtual Dependencies", [])]
        variables = hf["Variables"][()]
        if "Properties" in hf:
            result["properties"] = hf["Properties"][()]
    logger.info("  Infiltration model type: %s (%d classes)", result["model_type"], len(variables))

    cell_points = mesh.cell_centers[~mesh.cell_is_boundary]
    face_points = (mesh.face_points[mesh.faces_fp_idx[:, 0]]
                   + mesh.face_points[mesh.faces_fp_idx[:, 1]]) / 2
    layers = [inf_hdf.parent / d.replace("\\", "/").lstrip("./") for d in deps] or [inf_hdf]
    cell_cls = np.zeros(len(cell_points), dtype=np.int64)
    face_cls = np.zeros(len(face_points), dtype=np.int64)
    rasters = [_ClassRaster(layer) for layer in layers]

    def classes_at(points, cell_size=None):
        cls = np.zeros(len(points), dtype=np.int64)
        for r in rasters:
            cls = cls * len(r.ids) + r.classes_at(points, cell_size)
        return cls

    cell_cls = classes_at(cell_points)
    face_cls = classes_at(face_points)
    if cell_cls.max(initial=0) >= len(variables) or face_cls.max(initial=0) >= len(variables):
        raise ValueError(f"{inf_hdf}: classes from {', '.join(map(str, layers))} "
                         f"don't fit its {len(variables)} Variables rows")

    names = [n for n in variables.dtype.names if n != "Name"]
    result["face_classifications"] = face_cls.astype(np.int32)
    if composite:
        var_names = [n.decode("utf-8", "replace").strip() for n in variables["Name"]]
        nd = var_names.index("NoData") if "NoData" in var_names else 0
        table = {n: variables[n].astype(np.float32) for n in names}
        cls, avg = _composite_cells(mesh, {"inf": (classes_at, table, nd, np.arange(len(variables)))})["inf"]
        result["cell_classifications"] = cls
        result["variables"] = avg
    else:
        result["cell_classifications"] = cell_cls.astype(np.int32)
        result["variables"] = {n: variables[n][cell_cls].astype(np.float32) for n in names}
    result["available"] = True
    return result


# ============================================================================
# WORKFLOW B: READ MESH FROM EXISTING g01.hdf  (RASMapper-created mesh)
# ============================================================================

def read_mesh_from_g01hdf(g01_hdf_path: Path) -> tuple:
    """Read mesh topology from an existing g01.hdf created by RASMapper.

    RASMapper writes the Voronoi mesh topology (cells, faces, face points,
    connectivity) into g01.hdf when the user generates the 2D mesh.  The
    hydraulic tables (volume-elevation curves, face area-elevation) are NOT
    present yet – those are only added when "Compute Geometry" is clicked.

    Returns
    -------
    (mesh: MeshData, meta: dict)
        meta keys: area_name, default_mannings, dx, dy, title, projection_wkt
    """
    logger.info("Reading mesh topology from existing g01.hdf: %s", g01_hdf_path)

    with h5py.File(str(g01_hdf_path), 'r') as hf:
        # ---------- area name from Attributes ----------
        attrs_data = hf['Geometry/2D Flow Areas/Attributes'][()]
        area_name = attrs_data['Name'][0].decode('utf-8', errors='replace').strip('\x00').strip()
        default_mannings = float(attrs_data['Mann'][0])
        dx = float(attrs_data['Spacing dx'][0])
        dy = float(attrs_data['Spacing dy'][0])

        title = hf['Geometry'].attrs.get('Title', b'').decode('utf-8', errors='replace')
        projection_wkt = hf.attrs.get('Projection', b'').decode('utf-8', errors='replace')

        area_path = f'Geometry/2D Flow Areas/{area_name}'

        # ---------- seed points (interior cells) ----------
        cell_pts_input = hf['Geometry/2D Flow Areas/Cell Points'][()]       # (n_interior, 2)
        perimeter = hf[f'{area_path}/Perimeter'][()]                         # (n_perim_pts, 2)

        # ---------- full cell set (interior + ghost) ----------
        cell_centers = hf[f'{area_path}/Cells Center Coordinate'][()]       # (n_cells, 2)
        cells_face_info = hf[f'{area_path}/Cells Face and Orientation Info'][()]  # (n_cells, 2)
        cells_fp_idx_raw = hf[f'{area_path}/Cells FacePoint Indexes'][()]   # (n_cells, max_fp)

        # ---------- faces ----------
        faces_ci = hf[f'{area_path}/Faces Cell Indexes'][()]                # (n_faces, 2)
        faces_fp = hf[f'{area_path}/Faces FacePoint Indexes'][()]           # (n_faces, 2)
        faces_normal = hf[f'{area_path}/Faces NormalUnitVector and Length'][()]  # (n_faces, 3)

        # ---------- face points ----------
        fp_coords = hf[f'{area_path}/FacePoints Coordinate'][()]            # (n_fp, 2)
        fp_is_perim_raw = hf[f'{area_path}/FacePoints Is Perimeter'][()]    # (n_fp,) -1=perim, 0=interior
        if f'{area_path}/Faces Perimeter Info' in hf:
            perim_info = hf[f'{area_path}/Faces Perimeter Info'][()]
            perim_values = hf[f'{area_path}/Faces Perimeter Values'][()]
            faces_perimeter_pts = [perim_values[s:s + n] for s, n in perim_info]
        else:
            faces_perimeter_pts = [np.empty((0, 2))] * len(faces_ci)

    n_interior = len(cell_pts_input)
    n_cells = len(cell_centers)

    # Ghost cells are all cells whose index is >= n_interior (the seed-point count).
    # In standard RASMapper HDFs ghost cells also happen to have face_count == 1,
    # but in other formats (e.g. internal-tool HDFs) ghost cells can have more faces.
    cell_is_boundary = np.arange(n_cells) >= n_interior

    # Normalise fp_is_perimeter: True if perimeter (stored as -1 in HDF, 0 elsewhere)
    fp_is_perimeter = fp_is_perim_raw < 0  # bool array, True = perimeter

    # Fallback for internal-tool HDFs where FacePoints Is Perimeter is all zeros:
    # derive perimeter face points from topology — any FP belonging to a face
    # that borders a ghost (boundary) cell is a perimeter face point.
    if not fp_is_perimeter.any():
        ghost_cell_set = set(int(i) for i in np.where(cell_is_boundary)[0])
        perim_fp_set: set = set()
        for fi in range(len(faces_ci)):
            if int(faces_ci[fi, 0]) in ghost_cell_set or int(faces_ci[fi, 1]) in ghost_cell_set:
                perim_fp_set.add(int(faces_fp[fi, 0]))
                perim_fp_set.add(int(faces_fp[fi, 1]))
        fp_is_perimeter = np.array([i in perim_fp_set for i in range(len(fp_coords))],
                                   dtype=bool)

    mesh = MeshData()
    mesh.cell_centers = cell_centers.astype(np.float64)
    mesh.cell_points_input = cell_pts_input.astype(np.float64)
    mesh.face_points = fp_coords.astype(np.float64)
    mesh.perimeter = perimeter.astype(np.float64)
    mesh.cells_fp_idx = cells_fp_idx_raw.astype(np.int32)     # already padded with -1
    mesh.faces_cell_idx = faces_ci.astype(np.int32)
    mesh.faces_fp_idx = faces_fp.astype(np.int32)
    mesh.faces_normal_length = faces_normal.astype(np.float32)
    mesh.fp_is_perimeter = fp_is_perimeter.astype(np.int32)   # 1=perimeter, 0=interior
    mesh.cell_is_boundary = cell_is_boundary
    mesh.faces_perimeter_pts = faces_perimeter_pts

    meta = {
        "area_name": area_name,
        "default_mannings": default_mannings,
        "dx": dx,
        "dy": dy,
        "title": title,
        "projection_wkt": projection_wkt,
    }

    logger.info("  Mesh loaded: %d cells (%d interior, %d ghost), %d faces, %d face points",
                mesh.n_cells, mesh.n_interior, mesh.n_cells - mesh.n_interior,
                mesh.n_faces, mesh.n_fp)
    return mesh, meta


def compute_bc_external_faces(g01_hdf_path: Path, mesh: MeshData,
                               cell_dx: float = 500.0) -> np.ndarray:
    """Compute the Boundary Condition Lines / External Faces structured array.

    Algorithm:
      For each BC line polyline, HEC-RAS finds the contiguous chain of
      perimeter faces that spans from the face nearest the BC start endpoint
      to the face nearest the BC end endpoint (following perimeter adjacency).

      1. Build a perimeter-face adjacency graph (faces sharing a perimeter FP).
      2. For each BC line, find start-face (closest to BC line start) and
         end-face (closest to BC line end).
      3. BFS/walk along the perimeter from start-face to end-face in the
         direction that follows the BC line.
      4. Compute cumulative station = sum of face lengths along the chain.

    Returns a structured array with dtype:
      [('BC Line ID', i4), ('Face Index', i4),
       ('FP Start Index', i4), ('FP End Index', i4),
       ('Station Start', f4), ('Station End', f4)]
    """
    ext_dtype = np.dtype([
        ('BC Line ID', '<i4'), ('Face Index', '<i4'),
        ('FP Start Index', '<i4'), ('FP End Index', '<i4'),
        ('Station Start', '<f4'), ('Station End', '<f4'),
    ])

    with h5py.File(str(g01_hdf_path), 'r') as hf:
        if 'Geometry/Boundary Condition Lines/Polyline Points' not in hf:
            return np.array([], dtype=ext_dtype)
        poly_pts = hf['Geometry/Boundary Condition Lines/Polyline Points'][()]
        poly_info = hf['Geometry/Boundary Condition Lines/Polyline Info'][()]
        bc_attrs = hf['Geometry/Boundary Condition Lines/Attributes'][()]
        n_bc = len(poly_info)

    if n_bc == 0:
        return np.array([], dtype=ext_dtype)

    n_interior = mesh.n_interior
    faces_ci = mesh.faces_cell_idx
    faces_fp = mesh.faces_fp_idx
    fp_coords = mesh.face_points
    faces_nl = mesh.faces_normal_length    # (n_faces, 3) – column 2 is face length

    # --- Identify perimeter faces and their midpoints ---
    # Use cell_is_boundary directly to handle non-contiguous ghost cell indices
    # (internal-tool HDFs may have ghost cells at indices < n_interior)
    cell_is_bnd = mesh.cell_is_boundary
    perim_mask = cell_is_bnd[faces_ci[:, 0]] | cell_is_bnd[faces_ci[:, 1]]
    perim_face_indices = np.where(perim_mask)[0]
    midpoints = (fp_coords[faces_fp[perim_face_indices, 0]] +
                 fp_coords[faces_fp[perim_face_indices, 1]]) / 2.0

    # --- Build perimeter face adjacency (shared perimeter face point) ---
    # For each perimeter face point, list the perimeter faces that use it
    fp_to_perim_faces = {}  # fp_index -> list of perimeter-face-indices
    for pfi in perim_face_indices:
        for fp_idx in (int(faces_fp[pfi, 0]), int(faces_fp[pfi, 1])):
            fp_to_perim_faces.setdefault(fp_idx, []).append(int(pfi))

    # Adjacency: perim_face_adj[pfi] = list of adjacent perimeter faces
    perim_face_adj = {int(pfi): set() for pfi in perim_face_indices}
    for fp_idx, flist in fp_to_perim_faces.items():
        for i in range(len(flist)):
            for j in range(i + 1, len(flist)):
                perim_face_adj[flist[i]].add(flist[j])
                perim_face_adj[flist[j]].add(flist[i])

    all_rows = []
    for bc_id in range(n_bc):
        pt_start = int(poly_info[bc_id, 0])
        pt_count = int(poly_info[bc_id, 1])
        bc_pts = poly_pts[pt_start:pt_start + pt_count].astype(np.float64)
        bc_start = bc_pts[0]
        bc_end = bc_pts[-1]

        # --- Find start and end perimeter faces (closest to BC endpoints) ---
        def nearest_perim_face(target_pt):
            best_fi, best_dist = None, np.inf
            for idx, (pfi, mid) in enumerate(zip(perim_face_indices, midpoints)):
                d = np.linalg.norm(mid - target_pt)
                if d < best_dist:
                    best_dist, best_fi = d, int(pfi)
            return best_fi

        start_fi = nearest_perim_face(bc_start)
        end_fi = nearest_perim_face(bc_end)

        # --- Refined start/end using nearest PERIMETER FACE POINT (not midpoint) ---
        # When two adjacent faces share a face point near the BC endpoint,
        # HEC-RAS includes both (the face "before" and "after" that point).
        # Strategy: find nearest perimeter FP to BC start → find face containing it
        # that is earliest in the BC direction; similarly for BC end.

        fp_is_perim_arr = mesh.fp_is_perimeter   # 1=perimeter, 0=interior

        # Identify perimeter face points
        perim_fp_indices = np.where(fp_is_perim_arr > 0)[0]
        perim_fp_coords_arr = fp_coords[perim_fp_indices]

        def nearest_perim_fp(target_pt):
            dists = np.linalg.norm(perim_fp_coords_arr - target_pt, axis=1)
            return int(perim_fp_indices[np.argmin(dists)])

        def faces_containing_fp(fp_idx):
            """Perimeter faces that contain fp_idx."""
            return [int(pfi) for pfi in perim_face_indices
                    if fp_idx in (int(faces_fp[pfi, 0]), int(faces_fp[pfi, 1]))]

        start_fp_idx = nearest_perim_fp(bc_start)
        end_fp_idx = nearest_perim_fp(bc_end)

        # BC line direction for ordering candidates
        bc_dir = (bc_end - bc_start) / (np.linalg.norm(bc_end - bc_start) + 1e-9)

        def pick_start_face(fp_idx):
            """Face at start_fp that is the FIRST face in the BC chain.
            Rule: among faces containing fp_idx, exclude any whose OTHER face
            point has a NEGATIVE station (i.e. a spur/corner that goes backward
            past the BC start).  Among the remaining candidates pick the one
            whose other face point has the MINIMUM (earliest) positive station
            — that face is literally closest to the BC start along the perimeter.
            """
            cands = faces_containing_fp(fp_idx)
            if not cands:
                return nearest_perim_face(bc_start)
            def other_st(pfi):
                fp0, fp1 = faces_fp[pfi]
                other = fp1 if fp0 == fp_idx else fp0
                return float(np.dot(fp_coords[other] - bc_start, bc_dir))
            # exclude backward faces (other FP before BC start)
            fwd = [pfi for pfi in cands if other_st(pfi) >= 0]
            if not fwd:
                fwd = cands          # fallback: keep all
            return min(fwd, key=other_st)

        def pick_end_face(fp_idx):
            """Face at end_fp that is the LAST face in the BC chain.
            Rule: among faces containing fp_idx, pick the one whose OTHER face
            point has the MAXIMUM station — that is the face whose non-shared
            endpoint is furthest along the BC direction, meaning it is the
            'outermost' end of the chain (may extend slightly past the BC
            endpoint, which matches HEC-RAS behaviour)."""
            cands = faces_containing_fp(fp_idx)
            if not cands:
                return nearest_perim_face(bc_end)
            def other_st(pfi):
                fp0, fp1 = faces_fp[pfi]
                other = fp1 if fp0 == fp_idx else fp0
                return float(np.dot(fp_coords[other] - bc_start, bc_dir))
            return max(cands, key=other_st)

        start_fi = pick_start_face(start_fp_idx)
        end_fi = pick_end_face(end_fp_idx)

        if start_fi is None or end_fi is None:
            logger.warning("  BC line %d: could not find start/end faces", bc_id)
            continue

        chain = _bc_arc(start_fi, end_fi, perim_face_adj, bc_pts, fp_coords, faces_fp)
        if chain is None:
            logger.warning("  BC line %d: no path between start/end faces", bc_id)
            continue

        bc_name = bc_attrs['Name'][bc_id].decode('utf-8', errors='replace').strip('\x00')
        logger.info("  BC line %d (%s): %d faces in chain", bc_id, bc_name, len(chain))

        # --- Build External Faces with cumulative station (sum of face lengths) ---
        cum_station = 0.0
        # Orient each face so the chain is linked (end point of one face is
        # the start point of the next), as RASMapper writes it.
        first = [int(v) for v in faces_fp[chain[0]]]
        start_pt = first[0]
        if len(chain) > 1:
            not_shared = set(first) - set(int(v) for v in faces_fp[chain[1]])
            if len(not_shared) == 1:
                start_pt = not_shared.pop()
        for fi in chain:
            face_len = float(faces_nl[fi, 2])
            fa, fb = int(faces_fp[fi, 0]), int(faces_fp[fi, 1])
            fp1, fp2 = (fa, fb) if fa == start_pt else (fb, fa)
            start_pt = fp2
            all_rows.append((bc_id, fi, fp1, fp2,
                             np.float32(cum_station),
                             np.float32(cum_station + face_len)))
            cum_station += face_len

    if not all_rows:
        return np.array([], dtype=ext_dtype)
    return np.array(all_rows, dtype=ext_dtype)


def _write_infiltration_datasets(grp, inf_data: dict, n_interior: int,
                                  n_faces: int, gz: dict):
    """Write the infiltration datasets RASMapper writes: the cell and face
    classes, one dataset per variable of the layer, and its Properties.
    Uses _add_ds() so datasets already in the file are kept."""
    _add_ds(grp, "Cell Center Classifications", inf_data["cell_classifications"], gz)
    _add_ds(grp, "Face Center Classifications", inf_data["face_classifications"], gz)
    for name, values in inf_data["variables"].items():
        _add_ds(grp, name, values, gz)
    if inf_data.get("properties") is not None:
        _add_ds(grp, "Properties", inf_data["properties"], gz)


def update_geometry_hdf_with_tables(source_hdf: Path, output_path: Path,
                                     mesh: MeshData, cell_props: dict,
                                     face_props: dict, lc_data: dict,
                                     inf_data: dict, bc_ext_faces: np.ndarray,
                                     projection_wkt: str = ""):
    """Copy the RASMapper-created g01.hdf to output_path, then add the
    hydraulic tables (volume-elevation curves, face area-elevation profiles,
    Manning's n per cell, infiltration, percent-impervious, BC external faces)
    that the GUI computes when "Compute Geometry" is clicked.

    The mesh topology datasets already present in the source HDF are left
    intact; only the *missing* datasets / attributes are added.
    """
    import shutil
    logger.info("Updating geometry HDF with hydraulic tables: %s", output_path)

    if output_path != source_hdf:
        shutil.copy2(str(source_hdf), str(output_path))

    now_date = datetime.now().strftime("%d%b%Y %H:%M:%S").upper()
    now_time = datetime.now().strftime("%d%b%Y %H:%M:%S")

    n_cells = mesh.n_cells
    n_interior = mesh.n_interior
    n_faces = mesh.n_faces
    gz = {"compression": "gzip", "compression_opts": 1}

    with h5py.File(str(output_path), 'a') as hf:
        # --- Root attrs ---
        hf.attrs['File Type'] = np.bytes_('HEC-RAS Results')
        if projection_wkt:
            hf.attrs['Projection'] = np.bytes_(projection_wkt)

        # --- Geometry attrs ---
        geo = hf['Geometry']
        geo.attrs['Complete Geometry'] = np.bytes_('True')
        geo.attrs['Geometry Time'] = np.bytes_(now_time)

        # Detect area name
        attrs_data = hf['Geometry/2D Flow Areas/Attributes'][()]
        area_name = attrs_data['Name'][0].decode('utf-8', errors='replace').strip('\x00').strip()
        default_mannings = float(attrs_data['Mann'][0])
        area_path = f'Geometry/2D Flow Areas/{area_name}'
        ag = hf[area_path]

        # --- Perimeter 1 group attrs (added during Compute Geometry) ---
        ag.attrs['Cell Volume Tolerance'] = np.float32(0.01)
        ag.attrs['Cell Minimum Area Fraction'] = np.float32(0.01)
        ag.attrs['Composite LC'] = np.uint8(0)
        ag.attrs['Connection Profile Hash'] = np.bytes_('')
        ag.attrs['Face Area Conveyance Ratio'] = np.float32(0.02)
        ag.attrs['Face Area Elevation Tolerance'] = np.float32(0.01)
        ag.attrs['Face Profile Tolerance'] = np.float32(0.01)
        ag.attrs["Manning's n"] = np.float32(default_mannings)
        ag.attrs['Multiple Face Mann n'] = np.uint8(0)
        ag.attrs['Laminar Depth'] = np.float32(0.2)
        ag.attrs['Min Face Length Ratio'] = np.float32(np.nan)
        ag.attrs['Property Tables Last Computed'] = np.bytes_(now_time)
        ag.attrs['Property Tables LC Hash'] = np.bytes_('')
        ag.attrs['Terrain File Date'] = np.bytes_(now_date)
        ag.attrs['Terrain Filename'] = np.bytes_(r'.\Terrain\Terrain.hdf')
        ag.attrs['Land Cover Date Last Modified'] = np.bytes_(now_date)
        ag.attrs['Land Cover File Date'] = np.bytes_(now_date)
        ag.attrs['Land Cover Filename'] = np.bytes_(r'.\Land Classification\LandCover.hdf')
        ag.attrs['Land Cover Layername'] = np.bytes_('LandCover')
        ag.attrs['Infiltration Date Last Modified'] = np.bytes_(now_date)
        ag.attrs['Infiltration File Date'] = np.bytes_(now_date)
        ag.attrs['Infiltration Filename'] = np.bytes_(r'.\Land Classification\Infiltration.hdf')
        ag.attrs['Infiltration Layername'] = np.bytes_('Infiltration')
        ag.attrs['Percent Impervious Date Last Modified'] = np.bytes_(now_date)
        ag.attrs['Percent Impervious File Date'] = np.bytes_(now_date)
        ag.attrs['Percent Impervious Filename'] = np.bytes_(r'.\Land Classification\LandCover.hdf')
        ag.attrs['Percent Impervious Layername'] = np.bytes_('LandCover')

        # --- Per-cell hydraulic table datasets ---
        _add_ds(ag, 'Cells Center Manning\'s n',
                lc_data['mannings_n'].astype(np.float32), gz)
        _add_ds(ag, 'Cells Minimum Elevation',
                cell_props['min_elevations'], gz)
        _add_ds(ag, 'Cells Surface Area',
                cell_props['surface_areas'], gz)
        _add_ds(ag, 'Cells Volume Elevation Info',
                cell_props['vol_elev_info'], gz)
        _add_ds(ag, 'Cells Volume Elevation Values',
                cell_props['vol_elev_values'], gz)

        # --- Per-face hydraulic table datasets ---
        _add_ds(ag, 'Faces Minimum Elevation',
                face_props['min_elevations'], gz)
        _add_ds(ag, 'Faces Low Elevation Centroid',
                face_props['low_elev_centroid'], gz)
        _add_ds(ag, 'Faces Area Elevation Info',
                face_props['face_ae_info'], gz)
        _add_ds(ag, 'Faces Area Elevation Values',
                face_props['face_ae_values'], gz)

        # --- Infiltration subgroup ---
        # Only write when: the source HDF already has an Infiltration group
        # (already-computed project) OR an Infiltration.hdf was successfully
        # read (inf_data["available"] == True).  Do NOT create an Infiltration
        # group for projects that have no infiltration model (e.g. VA example)
        # because doing so changes simulation behaviour.
        source_has_inf = 'Infiltration' in ag
        if source_has_inf or inf_data.get("available", False):
            if 'Infiltration' not in ag:
                inf_grp = ag.create_group('Infiltration')
            else:
                inf_grp = ag['Infiltration']
            inf_grp.attrs['Infiltration Date Last Modified'] = np.bytes_(now_date)
            inf_grp.attrs['Infiltration File Date'] = np.bytes_(now_date)
            inf_grp.attrs['Infiltration Filename'] = np.bytes_(
                r'.\Land Classification\Infiltration.hdf')
            inf_grp.attrs['Infiltration Layername'] = np.bytes_('Infiltration')
            _write_infiltration_datasets(inf_grp, inf_data, n_interior, n_faces, gz)

        # --- Percent Impervious subgroup ---
        if 'Percent Impervious' not in ag:
            pi_grp = ag.create_group('Percent Impervious')
        else:
            pi_grp = ag['Percent Impervious']
        pi_grp.attrs['Percent Impervious Date Last Modified'] = np.bytes_(now_date)
        pi_grp.attrs['Percent Impervious File Date'] = np.bytes_(now_date)
        pi_grp.attrs['Percent Impervious Filename'] = np.bytes_(
            r'.\Land Classification\LandCover.hdf')
        pi_grp.attrs['Percent Impervious Layername'] = np.bytes_('LandCover')
        # Percent Impervious per cell (interior only)
        _add_ds(pi_grp, 'Percent Impervious',
                _percent_impervious(lc_data, n_cells, n_interior), gz)
        _add_ds(pi_grp, 'Cell Center Classifications',
                lc_data['classifications'][:n_interior].astype(np.int32), gz)
        _add_ds(pi_grp, 'Face Center Classifications',
                lc_data['face_classifications'].astype(np.int32), gz)

        # --- BC Lines External Faces ---
        bc_grp = hf['Geometry/Boundary Condition Lines']
        # Update Attributes: add 'Length' field if not present
        old_attrs = bc_grp['Attributes'][()]
        if 'Length' not in old_attrs.dtype.names:
            new_dt = np.dtype([('Name', 'S32'), ('SA-2D', 'S16'),
                               ('Type', 'S8'), ('Length', '<f4')])
            n_bc = len(old_attrs)
            new_data = np.zeros(n_bc, dtype=new_dt)
            for nm in ('Name', 'SA-2D', 'Type'):
                new_data[nm] = old_attrs[nm]
            # Compute lengths from External Faces stations
            for bc_id in range(n_bc):
                rows = bc_ext_faces[bc_ext_faces['BC Line ID'] == bc_id]
                if len(rows) > 0:
                    new_data['Length'][bc_id] = np.float32(
                        max(abs(rows['Station End']).max(), abs(rows['Station Start']).max()))
            del bc_grp['Attributes']
            bc_grp.create_dataset('Attributes', data=new_data, **gz)

        if len(bc_ext_faces) > 0:
            # Always overwrite External Faces — some source g01.hdf files
            # (e.g. from internal tools) already have an empty (0,) dataset
            # which _add_ds would silently skip.
            if 'External Faces' in bc_grp:
                del bc_grp['External Faces']
            bc_grp.create_dataset('External Faces', data=bc_ext_faces, **gz)

    logger.info("  Geometry HDF updated with hydraulic tables")



def _bc_arc(start_fi, end_fi, adj, bc_pts, fp_coords, faces_fp):
    """Perimeter faces from start_fi to end_fi that follow the BC polyline.

    The perimeter is a ring, so two arcs join the start and end faces. Take
    the one that lies along the polyline (smallest mean distance from points
    sampled on the polyline to the arc's face midpoints). Taking the shorter
    arc instead fails for a BC line that runs most of the way round the area,
    whose ends are close together.
    """
    def walk(first_step):
        path, prev, cur = [start_fi], start_fi, first_step
        while True:
            if cur == end_fi and not (cur == start_fi and len(path) == 1):
                return path + ([cur] if cur != start_fi else [])
            if cur == start_fi:
                return path          # went all the way round (start == end)
            path.append(cur)
            nxt = [n for n in adj[cur] if n != prev]
            if len(nxt) != 1:
                return None          # not a simple ring here
            prev, cur = cur, nxt[0]

    nbs = sorted(adj[start_fi])
    arcs = [[start_fi]] if start_fi == end_fi else []
    if len(nbs) == 2:
        for nb in nbs:
            a = walk(nb)
            if a:
                arcs.append(a)
    if not arcs:
        return _bfs_path(start_fi, end_fi, adj)

    seg = np.asarray(bc_pts, dtype=np.float64)
    lens = np.hypot(*np.diff(seg, axis=0).T)
    t = np.linspace(0, lens.sum(), 64)
    cum = np.r_[0, np.cumsum(lens)]
    samples = np.c_[np.interp(t, cum, seg[:, 0]), np.interp(t, cum, seg[:, 1])]

    def cost(arc):
        mids = (fp_coords[faces_fp[arc, 0]] + fp_coords[faces_fp[arc, 1]]) / 2
        return np.sqrt(((samples[:, None, :] - mids[None]) ** 2).sum(-1)).min(1).mean()

    return min(arcs, key=cost)


def _bfs_path(start_fi, end_fi, adj):
    from collections import deque
    visited = {start_fi: None}
    queue = deque([start_fi])
    while queue:
        cur = queue.popleft()
        if cur == end_fi:
            path = []
            while cur is not None:
                path.append(cur)
                cur = visited[cur]
            return path[::-1]
        for nb in sorted(adj[cur]):
            if nb not in visited:
                visited[nb] = cur
                queue.append(nb)
    return None

def _add_ds(grp, name: str, data, gz: dict):
    """Add dataset to group only if it does not already exist."""
    if name not in grp:
        grp.create_dataset(name, data=data, **gz)
    # If it exists, leave it (already computed by GUI or previous run)


def _percent_impervious(lc_data, n_cells, n_interior):
    """Percent impervious per interior cell: its cell average under Composite
    LC, otherwise its land-cover class's value."""
    if lc_data.get("percent_impervious") is not None:
        return lc_data["percent_impervious"][:n_interior].astype(np.float32)
    id_to_name = lc_data.get('id_to_name', {})
    name_to_imperv = lc_data.get('name_to_imperv', {})
    lc_class = lc_data.get('classifications', np.zeros(n_cells, dtype=np.int32))
    return np.array([name_to_imperv.get(id_to_name.get(int(lc_class[ci]), ''), 0.0)
                     for ci in range(n_interior)], dtype=np.float32)


# ============================================================================
# STEP 7: ASSEMBLE p01.tmp.hdf
# ============================================================================

def assemble_p01_tmp_hdf(g01_hdf_path: Path, p01_text_path: Path,
                          output_path: Path, projection_wkt: str = "",
                          project_name: str = "",
                          event_conditions_from: Path = None):
    """Assemble p01.tmp.hdf from generated g01.hdf + parsed .p01 text.

    Event Conditions are empty unless ``event_conditions_from`` names an HDF
    (a plan HDF such as a GUI-run p01.hdf, or a u01.hdf) whose
    "Event Conditions" group is copied in. Boundary conditions in it are then
    re-pointed at this geometry's faces (see remap_bc_faces).
    """
    logger.info("Assembling p01.tmp.hdf...")

    with h5py.File(str(g01_hdf_path), 'r') as fg:
        with h5py.File(str(output_path), 'w') as fd:
            # Copy root attrs
            for k, v in fg.attrs.items():
                fd.attrs[k] = v

            # Copy Geometry
            fg.copy('Geometry', fd)

            copied = False
            if event_conditions_from:
                with h5py.File(str(event_conditions_from), 'r') as fe:
                    if 'Event Conditions' in fe:
                        fe.copy('Event Conditions', fd)
                        copied = True
                        logger.info("  Event Conditions copied from %s", event_conditions_from)
                    else:
                        logger.warning("  %s has no Event Conditions group", event_conditions_from)
            if not copied:
                ec = fd.create_group('Event Conditions')
                ec_u = ec.create_group('Unsteady')
                ec_u.create_group('Boundary Conditions')
                ec_u.create_group('Initial Conditions')

            # Plan Data from .p01 text
            _write_plan_data(fd, p01_text_path, project_name=project_name)

    if copied:
        if _same_mesh(event_conditions_from, output_path):
            logger.info("  Event Conditions come from this mesh; BC face lists kept as they are")
        else:
            remap_bc_faces(output_path)
    logger.info("  p01.tmp.hdf assembled: %s", output_path)


def _same_mesh(hdf_a: Path, hdf_b: Path) -> bool:
    """True if both files hold the same 2D mesh (faces and face points)."""
    keys = ("Faces FacePoint Indexes", "FacePoints Coordinate")
    try:
        with h5py.File(str(hdf_a), 'r') as fa, h5py.File(str(hdf_b), 'r') as fb:
            ga, gb = fa['Geometry/2D Flow Areas'], fb['Geometry/2D Flow Areas']
            areas = [k for k in gb if isinstance(gb[k], h5py.Group)]
            for area in areas:
                if area not in ga:
                    return False
                for k in keys:
                    a, b = ga[area][k][()], gb[area][k][()]
                    if a.shape != b.shape or not np.allclose(a, b):
                        return False
            return bool(areas)
    except (KeyError, OSError):
        return False


def remap_bc_faces(plan_hdf: Path) -> int:
    """Point each 2D boundary condition in Event Conditions at this geometry.

    The solver applies a 2D BC to the faces listed in the BC dataset's
    "Face Indexes" / "Face Point Indexes" / "Face Fraction" attributes, not to
    Geometry/Boundary Condition Lines/External Faces. Those lists belong to
    the mesh they were computed on, so BC data copied from another run must be
    re-pointed or the solver reports "overlapping boundaries" or applies flow
    to the wrong faces. Faces are taken whole (fraction 1); RASMapper may use
    a partial fraction for the last face of a line.
    Returns the number of BC datasets updated.
    """
    n = 0
    with h5py.File(str(plan_hdf), 'a') as hf:
        bcl = hf.get('Geometry/Boundary Condition Lines')
        bcs = hf.get('Event Conditions/Unsteady/Boundary Conditions')
        if bcl is None or bcs is None or 'External Faces' not in bcl:
            return 0
        names = [x.decode('utf-8', 'replace').strip('\x00').strip()
                 for x in bcl['Attributes'][()]['Name']]
        ext = bcl['External Faces'][()]

        def visit(name, obj):
            nonlocal n
            if not isinstance(obj, h5py.Dataset) or 'BC Line' not in obj.attrs:
                return
            line = obj.attrs['BC Line']
            line = line.decode('utf-8', 'replace') if isinstance(line, bytes) else str(line)
            line = line.strip('\x00').strip()
            if line not in names:
                logger.warning("  BC %r: no BC line of that name in the geometry", line)
                return
            rows = ext[ext['BC Line ID'] == names.index(line)]
            if len(rows) == 0:
                logger.warning("  BC %r: no perimeter faces assigned", line)
                return
            obj.attrs['Face Indexes'] = rows['Face Index'].astype(np.int32)
            obj.attrs['Face Point Indexes'] = np.r_[rows['FP Start Index'][:1],
                                                    rows['FP End Index']].astype(np.int32)
            obj.attrs['Face Fraction'] = np.ones(len(rows), dtype=np.float32)
            n += 1
        bcs.visititems(visit)
    logger.info("  Re-pointed %d boundary condition(s) at this geometry's faces", n)
    return n


def _write_plan_data(hf, p01_path, project_name: str = ""):
    """Parse .p01 text and write Plan Data group."""
    raw = {}
    with open(p01_path, 'r') as f:
        for line in f:
            line = line.rstrip('\n\r')
            if '=' in line:
                key, val = line.split('=', 1)
                raw[key.strip()] = val.strip()

    # Use project_name for filenames; fall back to p01_path stem if not provided
    pname = project_name or Path(p01_path).stem.split('.')[0]

    pd = hf.create_group('Plan Data')
    pi = pd.create_group('Plan Information')

    pi.attrs['Plan Title'] = np.bytes_(raw.get('Plan Title', ''))
    pi.attrs['Plan Name'] = np.bytes_(raw.get('Plan Title', ''))
    pi.attrs['Plan ShortID'] = np.bytes_(raw.get('Short Identifier', '01').strip())
    pi.attrs['Computation Time Step Base'] = np.bytes_(raw.get('Computation Interval', '1MIN'))
    pi.attrs['Base Output Interval'] = np.bytes_(raw.get('Output Interval', '1HOUR'))
    pi.attrs['Geometry Filename'] = np.bytes_(f"{pname}.{raw.get('Geom File', 'g01')}")
    pi.attrs['Flow Filename'] = np.bytes_(f"{pname}.{raw.get('Flow File', 'u01')}")
    pi.attrs['Plan Filename'] = np.bytes_(f"{pname}.{raw.get('Plan Title', 'p01')}")
    pi.attrs['Project Title'] = np.bytes_(raw.get('Plan Title', 'Project'))
    pi.attrs['Project Filename'] = np.bytes_('project.prj')

    # Parse simulation dates
    sim_date = raw.get('Simulation Date', '')
    if sim_date:
        parts = sim_date.split(',')
        if len(parts) >= 4:
            st = f"{parts[0]} {parts[1][:2]}:{parts[1][2:]}:00"
            et = f"{parts[2]} {parts[3][:2]}:{parts[3][2:]}:00"
            pi.attrs['Simulation Start Time'] = np.bytes_(st)
            pi.attrs['Simulation End Time'] = np.bytes_(et)
            pi.attrs['Time Window'] = np.bytes_(f"{st} to {et}")

    # Plan Parameters — all with exact HDF dtypes
    pp = pd.create_group('Plan Parameters')
    _write_plan_parameters(pp, raw)


def _write_plan_parameters(pp, raw):
    """Write all 58 Plan Parameters with exact HDF dtypes."""
    def f32(v): return np.float32(float(v or '0'))
    def i32(v): return np.int32(int(v or '0'))
    def b(v): return np.bytes_(v)
    def af32(v): return np.array(v, dtype=np.float32)
    def ai32(v): return np.array(v, dtype=np.int32)
    def au8(v): return np.array(v, dtype=np.uint8)
    def ab(v): return np.array([np.bytes_(x) for x in v])

    # 1D params
    pp.attrs['1D Methodology'] = b(raw.get('UNET 1D Methodology', 'Finite Difference'))
    pp.attrs['1D Theta'] = f32(raw.get('UNET Theta', '1'))
    pp.attrs['1D Theta Warmup'] = f32(raw.get('UNET Theta Warmup', '1'))
    pp.attrs['1D Water Surface Elevation Tolerance'] = f32(raw.get('UNET ZTol', '0.02'))
    pp.attrs['1D Storage Area Elevation Tolerance'] = f32(raw.get('UNET ZSATol', '0.02'))
    pp.attrs['1D Flow Tolerance'] = np.float32('nan')
    pp.attrs['1D Maximum Iterations'] = i32(raw.get('UNET MxIter', '20'))
    pp.attrs['1D Maximum Iterations Without Improvement'] = i32(raw.get('UNET Max Iter WO Improvement', '0'))
    pp.attrs['1D Maximum Water Surface Error To Abort'] = f32(raw.get('UNET DZMax Abort', '100'))
    pp.attrs['1D Cores'] = i32(raw.get('UNET D1 Cores', '0'))

    # 1D-2D coupling
    pp.attrs['1D-2D Maximum Iterations'] = i32(raw.get('UNET D1D2 MaxIter', '0'))
    pp.attrs['1D-2D Water Surface Tolerance'] = f32(raw.get('UNET D1D2 ZTol', '0.01'))
    pp.attrs['1D-2D Flow Tolerance'] = f32(raw.get('UNET D1D2 QTol', '0.1'))
    pp.attrs['1D-2D Minimum Flow Tolerance'] = f32(raw.get('UNET D1D2 MinQTol', '1'))
    pp.attrs['1D-2D Number of Warmup Steps'] = i32(0)
    pp.attrs['1D-2D Warmup Time Step (hours)'] = f32(raw.get('UNET DtIC', '0'))
    pp.attrs['1D-2D Maximum Number of Time Slices'] = i32(raw.get('UNET MaxCRTS', '20'))
    pp.attrs['1D-2D Minimum Time Step for Slicing(hours)'] = f32(raw.get('UNET DtMin', '0'))
    pp.attrs['1D-2D Weir Flow Submergence Decay Exponent'] = f32(raw.get('UNET WFX', '1'))
    pp.attrs['1D-2D Gate Flow Submergence Decay Exponent'] = f32(raw.get('UNET SFX', '1'))
    pp.attrs['1D-2D IS Stablity Factor'] = f32(raw.get('UNET SFStab', '1'))
    pp.attrs['1D-2D LS Stablity Factor'] = f32(raw.get('UNET WFStab', '2'))

    # 2D params (arrays)
    eq_map = {'0': 'Diffusion Wave', '1': 'Shallow Water Equations'}
    pp.attrs['2D Equation Set'] = ab([eq_map.get(raw.get('UNET D2 Equation', '0'), 'Diffusion Wave')])
    pp.attrs['2D Theta'] = af32([float(raw.get('UNET D2 Theta', '1') or '1')])
    pp.attrs['2D Theta Warmup'] = af32([float(raw.get('UNET D2 Theta Warmup', '1') or '1')])
    pp.attrs['2D Water Surface Tolerance'] = af32([float(raw.get('UNET D2 Z Tol', '0.01') or '0.01')])
    pp.attrs['2D Volume Tolerance'] = af32([float(raw.get('UNET D2 Volume Tol', '0.01') or '0.01')])
    pp.attrs['2D Maximum Iterations'] = ai32([int(raw.get('UNET D2 Max Iterations', '20') or '20')])
    pp.attrs['2D Number of Time Slices'] = ai32([int(raw.get('UNET D2 TimeSlices', '1') or '1')])
    pp.attrs['2D Cores (per mesh)'] = ai32([int(raw.get('UNET D2 Cores', '0') or '0')])
    pp.attrs['2D Coriolis'] = b('True' if raw.get('UNET D2 Coriolis', '0') != '0' else 'False')
    pp.attrs['2D Latitude for Coriolis'] = af32([3.4028235e+38])
    pp.attrs['2D Turbulence Formulation'] = ab([raw.get('UNET D2 Turbulence Formulation', 'None')])
    pp.attrs['2D Smagorinsky Mixing Coefficient'] = af32([float(raw.get('UNET D2 Smagorinsky Mixing', '0') or '0')])
    pp.attrs['2D Longitudinal Mixing Coefficient'] = af32([float(raw.get('UNET D2 Eddy Viscosity', '0') or '0')])
    pp.attrs['2D Transverse Mixing Coefficient'] = af32([float(raw.get('UNET D2 Transverse Eddy Viscosity', '0') or '0')])
    pp.attrs['2D Initial Conditions Ramp Up Time (hrs)'] = af32([float(raw.get('UNET D2 TotalICTime', '0') or '0')])
    pp.attrs['2D Boundary Condition Ramp Up Fraction'] = af32([float(raw.get('UNET D2 RampUpFraction', '0.1') or '0.1')])
    pp.attrs['2D Boundary Condition Volume Check'] = ab(['True' if raw.get('UNET D2 BCVolumeCheck', '0') != '0' else 'False'])
    pp.attrs['2D Matrix Solver'] = ab([raw.get('UNET D2 SolverType', 'PARDISO (Direct)').split('(')[0].strip()])
    pp.attrs['2D Advanced Convergence'] = au8([int(raw.get('UNET D2 Advanced Convergence', '0') or '0')])
    pp.attrs['2D WS Max Tolerance'] = af32([float(raw.get('UNET D2 WS Max Tol', '0') or '0')])
    pp.attrs['2D WS RMS Tolerance'] = af32([float(raw.get('UNET D2 WS RMS Tol', '0') or '0')])
    pp.attrs['2D WS Stalling Tolerance'] = af32([float(raw.get('UNET D2 WS Stall Tol', '0') or '0')])
    pp.attrs['2D Only'] = b('True')
    pp.attrs['2D Names'] = ab(['Perimeter 1'])

    # Other
    pp.attrs['Gravity'] = f32(raw.get('UNET Gravity', '32.174'))
    pp.attrs['Friction Slope Average Method (XS)'] = b('Average Friction Slope')
    pp.attrs['Friction Slope Average Method (BR)'] = b('Average Conveyance')
    pp.attrs['Pardiso Solver'] = b('True' if raw.get('UNET Pardiso', '0') != '0' else 'False')
    pp.attrs['HDF Write Warmup'] = b('True' if raw.get('HDF Write Warmup', '0') != '0' else 'False')
    pp.attrs['HDF Write Time Slices'] = b('True' if raw.get('HDF Write Time Slices', '0') != '0' else 'False')
    pp.attrs['HDF Flush Buffer'] = b('True' if raw.get('HDF Flush', '0') != '0' else 'False')
    pp.attrs['HDF Compression'] = i32(raw.get('HDF Compression', '1'))
    pp.attrs['HDF Chunk Size'] = f32(raw.get('HDF Chunk Size', '1'))
    pp.attrs['HDF Spatial Parts'] = i32(raw.get('HDF Spatial Parts', '1'))
    pp.attrs['HDF Use Max Rows'] = i32(raw.get('HDF Use Max Rows', '0'))
    pp.attrs['HDF Fixed Rows'] = i32(raw.get('HDF Fixed Rows', '1'))


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description="HEC-RAS 2D Preprocessor for Linux")
    parser.add_argument("project_dir", type=Path, help="Path to HEC-RAS project directory")
    parser.add_argument("--project-name", default=None, help="Project name (auto-detected from .prj)")
    parser.add_argument("--plan", default="01", help="Plan number (default: 01)")
    parser.add_argument("--output-dir", type=Path, default=None, help="Output directory (default: project_dir)")
    parser.add_argument("--event-conditions", type=Path, default=None, metavar="HDF",
                        help="Workflow B: copy boundary/initial conditions from this HDF "
                             "(e.g. a GUI-run .p01.hdf) instead of leaving them to the DSS file")
    args = parser.parse_args()

    project_dir = args.project_dir
    output_dir = args.output_dir or project_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    # Auto-detect project name: prefer .g01 (unambiguous) over .prj
    # because some projects have auxiliary .prj (coordinate system) files
    # that do NOT share the project name (e.g. "Polygon Layer.prj").
    name = args.project_name
    if not name:
        g01_files = list(project_dir.glob("*.g01"))
        if g01_files:
            name = g01_files[0].stem
        else:
            # Fall back to HEC-RAS .prj: a real HEC-RAS .prj file has a
            # matching .g01 or .g01.hdf.  Skip any that don't.
            prj_files = list(project_dir.glob("*.prj"))
            for pf in prj_files:
                stem = pf.stem
                if (project_dir / f"{stem}.g01").exists() or \
                   (project_dir / f"{stem}.g01.hdf").exists():
                    name = stem
                    break
            if not name and prj_files:
                name = prj_files[0].stem   # last resort
        if not name:
            logger.error("Cannot detect project name. Use --project-name")
            sys.exit(1)

    plan = args.plan.zfill(2)
    logger.info("Project: %s, Plan: %s", name, plan)

    # Find files
    g01_path = project_dir / f"{name}.g01"
    p01_path = project_dir / f"{name}.p{plan}"
    terrain_tif = _find_terrain_tif(project_dir)
    geom_hdf = project_dir / f"{name}.g{plan}.hdf"
    lc_hdf = _find_layer(project_dir, geom_hdf, "Land Cover Filename", "LandCover.hdf")
    lc_tif = lc_hdf.with_suffix(".tif")
    inf_hdf = _find_layer(project_dir, geom_hdf, "Infiltration Filename", "Infiltration.hdf")

    # Check required files
    for f, desc in [(g01_path, ".g01"), (p01_path, f".p{plan}"), (terrain_tif, "Terrain TIF")]:
        if not f or not f.exists():
            logger.error("Required file not found: %s (%s)", f, desc)
            sys.exit(1)

    # Read projection from terrain
    ds = gdal.Open(str(terrain_tif))
    projection_wkt = ds.GetProjection() if ds else ""
    ds = None

    # -----------------------------------------------------------------------
    # Detect workflow:
    #   A) g01.hdf already has hydraulic tables AND p01.hdf (or p01.tmp.hdf)
    #      already exists → strip Results, copy files, skip Python preprocessing
    #   B) g01.hdf exists but no hydraulic tables (RASMapper created the mesh)
    #      → read mesh from HDF, compute hydraulic tables, update HDF in-place
    # -----------------------------------------------------------------------
    src_g01_hdf = project_dir / f"{name}.g{plan}.hdf"
    out_g01_hdf = output_dir / f"{name}.g{plan}.hdf"

    src_p01_hdf = project_dir / f"{name}.p{plan}.hdf"
    src_p01_tmp = project_dir / f"{name}.p{plan}.tmp.hdf"
    p01_tmp_path = output_dir / f"{name}.p{plan}.tmp.hdf"

    if src_g01_hdf.exists() and _g01hdf_has_tables(src_g01_hdf) and \
            (src_p01_hdf.exists() or src_p01_tmp.exists()):
        # ===================================================================
        # WORKFLOW A: GUI already ran geometry compute — just strip results
        # and hand off to RasGeomPreprocess + RasUnsteady.
        # ===================================================================
        logger.info("=== WORKFLOW A: Fully-computed g01.hdf + existing p01 found ===")
        logger.info("    Copying geometry HDF and stripping results from plan HDF …")

        # Copy g01.hdf to output (no-op when same path)
        if src_g01_hdf.resolve() != out_g01_hdf.resolve():
            shutil.copy2(str(src_g01_hdf), str(out_g01_hdf))
            logger.info("    Copied: %s", out_g01_hdf)

        # Produce p01.tmp.hdf: prefer pre-existing tmp.hdf; else strip p01.hdf
        if src_p01_tmp.exists():
            if src_p01_tmp.resolve() != p01_tmp_path.resolve():
                shutil.copy2(str(src_p01_tmp), str(p01_tmp_path))
                logger.info("    Copied: %s", p01_tmp_path)
        else:
            logger.info("    Stripping results from %s …", src_p01_hdf.name)
            _strip_results(src_p01_hdf, p01_tmp_path)
            logger.info("    Written: %s", p01_tmp_path)

        logger.info("=== DONE (Workflow A) ===")
        logger.info("  Geometry: %s", out_g01_hdf)
        logger.info("  Plan:     %s", p01_tmp_path)
        logger.info("")
        logger.info("Next steps:")
        logger.info("  1. RasGeomPreprocess %s.p%s.tmp.hdf x%s", name, plan, plan)
        logger.info("  2. RasUnsteady %s.p%s.tmp.hdf x%s", name, plan, plan)
        return

    if not src_g01_hdf.exists():
        logger.error("No %s: Workflow B needs the 2D mesh in the geometry HDF. "
                     "Build it in RASMapper, or with a mesh generator that writes one.",
                     src_g01_hdf.name)
        sys.exit(1)

    logger.info("=== WORKFLOW B: RASMapper g01.hdf found – computing hydraulic tables ===")

    # Step A1: Read mesh topology from existing HDF
    mesh, meta = read_mesh_from_g01hdf(src_g01_hdf)
    default_mannings = meta["default_mannings"]
    if not projection_wkt:
        projection_wkt = meta.get("projection_wkt", "")

    settings = table_settings(src_g01_hdf, meta["area_name"])

    # Step A2: Sample land cover
    lc_data = _default_land_cover(mesh, default_mannings)
    if lc_hdf.exists() and lc_tif and lc_tif.exists():
        lc_data = sample_land_cover(mesh, str(lc_hdf), str(lc_tif), default_mannings,
                                    composite=settings["composite_lc"])

    # Step A3: Face tables, then cell tables (a cell's minimum depends on its faces)
    terrain_path = _find_terrain_layer(project_dir, src_g01_hdf) or terrain_tif
    logger.info("Terrain: %s", terrain_path)
    terrain = property_tables.TerrainLayer(str(terrain_path))
    face_props = compute_face_properties(mesh, terrain, lc_data["face_mannings_n"], settings)
    cell_props = compute_cell_properties(mesh, terrain, face_props["min_elevations"], settings)

    # Step A5: Infiltration
    inf_data = sample_infiltration(mesh, str(inf_hdf) if inf_hdf.exists() else "",
                                   composite=settings["composite_lc"])

    # Step A6: Compute BC External Faces
    bc_ext_faces = compute_bc_external_faces(src_g01_hdf, mesh, meta.get("dx", 500.0))

    # Step A7: Copy source HDF to output_dir and add tables
    update_geometry_hdf_with_tables(src_g01_hdf, out_g01_hdf, mesh,
                                     cell_props, face_props, lc_data,
                                     inf_data, bc_ext_faces, projection_wkt)

    # Step 7: Assemble p01.tmp.hdf
    assemble_p01_tmp_hdf(out_g01_hdf, p01_path, p01_tmp_path, projection_wkt,
                         project_name=name, event_conditions_from=args.event_conditions)

    # Step 8: Generate .b01 runtime boundary conditions file
    u01_path = project_dir / f"{name}.u{plan}"
    b01_path = output_dir / f"{name}.b{plan}"
    generate_b01(p01_path, u01_path, b01_path, name)

    # Step 9: Copy .u01 and generate .x01 to output dir
    if u01_path.exists() and output_dir != project_dir:
        shutil.copy2(str(u01_path), str(output_dir / u01_path.name))
    x01_path = output_dir / f"{name}.x{plan}"
    if not x01_path.exists():
        _generate_x01(x01_path, out_g01_hdf)

    logger.info("=== DONE ===")
    logger.info("  Geometry: %s", out_g01_hdf)
    logger.info("  Plan:     %s", p01_tmp_path)
    logger.info("")
    logger.info("Next steps:")
    logger.info("  1. RasGeomPreprocess %s.p%s.tmp.hdf x%s", name, plan, plan)
    logger.info("  2. RasUnsteady %s.p%s.tmp.hdf x%s", name, plan, plan)


def generate_b01(p01_path: Path, u01_path: Path, b01_path: Path,
                  project_name: str):
    """Generate .b01 runtime boundary conditions file from .p01 and .u01.

    The .b01 is a text file (CRLF line endings) required by RasUnsteady.
    Format derived from HEC-RAS 6.6 GUI-generated .b01 files.
    """
    # Parse .p01
    raw = {}
    with open(p01_path, 'r') as f:
        for line in f:
            line = line.rstrip('\n\r')
            if '=' in line:
                key, val = line.split('=', 1)
                raw[key.strip()] = val.strip()

    plan_title = raw.get('Plan Title', project_name)
    short_id = raw.get('Short Identifier', plan_title).strip().ljust(64)

    # Parse simulation dates: "02JAN2000,0000,02JAN2000,2400"
    # HEC-RAS .b01 uses previous-day-at-2400 for 0000 start times
    sim_date = raw.get('Simulation Date', '')
    parts = sim_date.split(',') if sim_date else []
    if len(parts) >= 4:
        start_date_raw, start_time_raw = parts[0].strip(), parts[1].strip()
        end_date_raw, end_time_raw = parts[2].strip(), parts[3].strip()

        # Convert 0000 → previous day at 2400
        if start_time_raw == '0000':
            from datetime import datetime as dt, timedelta
            d = dt.strptime(start_date_raw.upper(), '%d%b%Y')
            d -= timedelta(days=1)
            start_b01 = d.strftime('%d%b%Y').lstrip('0') + ' 2400'
            # Fix: strftime may give '1Jan2000' or '01Jan2000' — ensure 2-digit day
            day_str = d.strftime('%d')
            mon_str = d.strftime('%b')
            yr_str = d.strftime('%Y')
            start_b01 = f"{day_str}{mon_str}{yr_str} 2400"
        else:
            start_b01 = f"{start_date_raw} {start_time_raw}"

        end_b01 = f"{end_date_raw.title()} {end_time_raw}"
    else:
        start_b01 = '01Jan2000 2400'
        end_b01 = '02Jan2000 2400'

    # Count actual BCs from .u01 (excluding empty pairs for FA)
    n_bc = 0
    if u01_path.exists():
        with open(u01_path, 'r') as f:
            for line in f:
                if line.startswith('Boundary Location='):
                    # Skip empty boundary location (pair for FA/precipitation)
                    fields = line.split(',')
                    bc_name = fields[7].strip() if len(fields) > 7 else ''
                    if bc_name:
                        n_bc += 1
    n_bc = max(n_bc, 2)  # minimum 2 (upstream + downstream)

    # Parse tolerances from .p01
    dz_tol = raw.get('UNET ZTol', '0.01').strip() or '0.01'
    dzsa_tol = raw.get('UNET ZSATol', '0.01').strip() or '0.01'
    warmup_steps = raw.get('UNET DtIC', '0').strip() or '0'
    # GUI uses 20 warmup steps by default for internal tool models
    warmup_int = int(float(warmup_steps)) if warmup_steps != '0' else 20

    dss_path = f"{project_name}.dss"

    # Build .b01 content (CRLF line endings)
    lines = [
        'HEC-RAS 6.6 September 2024',
        '       1       1       0       0',
        '       0       0',
        '       F',
        '       1',
        'Initial Conditions Flow Information',
        '       1        Initial Profile',
        ' 3.1E+38       2',
        'Flow and Seasonal Roughness Flag (plan)',
        '       F       F',
        '       2       0       0',
        '       3       0    .001',
        ' 3.1E+38',
        '       1       1       0       0',
        '       0       0',
        '       F',
        'Project Title, Plan Title and Plan ShortID',
        project_name,
        plan_title,
        short_id,
        'Job Control Information',
        '  Computation Interval  = %s' % raw.get('Computation Interval', '1MIN'),
        '  Warmup Interval       =  0 ',
        '  Instantaneous Profile = %s' % raw.get('Instantaneous Interval', '1HOUR'),
        '  Hydrograph Interval   = %s' % raw.get('Output Interval', '1HOUR'),
        '  Theta Simulation      =        1',
        '  Theta Warmup          =        1',
        '  Friction Slope Method =        2',
        '  Maximum Iterations    =       20',
        '  Max Iter WOImprovement=        0',
        '  Number Warmup Steps   =       %d' % warmup_int,
        '  Abort DZ Tolerance    =      100',
        '  DZ Tolerance          =      %s' % dz_tol.lstrip('0'),
        '  DZSA Tolerance        =      %s' % dzsa_tol.lstrip('0'),
        '  DQ Tolerance          =  3.1E+38',
        '  Weir Flow Stability   =        2',
        '  Spillway Stability    =        1',
        '  Write Restart File    =        F       F',
        '  Echo Input TS         =        F',
        '  Echo Parameters       =        F',
        '  Echo Output TS        =        F',
        '  DSS Message Level     =        4',
        '  Write HDF5 File       =        T',
        '  Write DSS File        =        T',
        dss_path,
        'Computational Time Window',
        '  Start Date/Time       = %s' % start_b01,
        '  End Date/Time         = %s' % end_b01,
        'Initial Conditions (use restart file?)',
        '       F',
        'Log File Information',
        '       F       0       0',
        'Computation Level Output',
        '       F-3.4E+38 3.4E+38       F       F       F       F       F       F       F       F       F   1HOUR',
        'Mixed Flow - Acceleration term reduction based on Froude number',
        '       F       4      .8       1',
        'Number of Gate Groups and Internal Boundaries with Gates',
        '       0       0       0       0',
        'Breach Data',
        '       0',
        'Hydrograph Data',
        '      %d' % n_bc,
        '       F       F       T       F       F               F',
        'Upstream Flow Hydrograph - River: Fake River  Reach: Fake Reach  RS: 100',
        '       2',
        '       0     100    8760     100',
        ' 3.4E+38',
        '       F       F       F       T       F',
        'Downstream Normal Depth',
        '    .001',
        'Internal Observed Stage/Flow Boundaries',
        '       0',
        'Ground Water Interflows',
        '       0',
        'Old River Diversions',
        '       F',
        'Lateral Inflows, Ungaged Lateral Inflows, Outlet TS, and Observed DSS',
        '       0       0       0       0',
        'Stage and Flow Boundary and Ungaged Areas',
        '       0       0',
        'Time Slicing Parameters',
        '       F',
        'HYDROGRAPH LOCATIONS',
        ' 0 ',
        'Rules (number of rule sets, number of lookbacks, number of tables)',
        '       0       0       0',
        'Extra Commands',
        '       0',
    ]

    with open(b01_path, 'wb') as f:
        f.write('\r\n'.join(lines).encode('latin-1'))

    logger.info("  Generated .b01: %s (%d BCs, time %s to %s)",
                b01_path.name, n_bc, start_b01, end_b01)


def _generate_x01(x01_path: Path, g01_hdf_path: Path):
    """Generate a minimal .x01 (cross-section index) file.

    Format derived from HEC-RAS 6.6 GUI-generated .x01 files.
    Uses BEC template with substituted area name and cell count.
    """
    with h5py.File(str(g01_hdf_path), 'r') as f:
        attrs = f['Geometry/2D Flow Areas/Attributes'][()]
        area_name = attrs['Name'][0].decode().strip('\x00').strip()
        # The storage-area record holds the number of perimeter points
        # (GUI-written .x01 files match "Polygon Points" in every example).
        n_perim = len(f['Geometry/2D Flow Areas/Polygon Points'])

    # Pad area name to 16 chars (HEC-RAS field width)
    area_padded = area_name.ljust(16)[:16]

    lines = [
        'HEC-RAS 6.6 September 2024',
        'Section - Arrays Sizes',
        '01',
        '       1       1       0       2       8       F',
        '       3       3       0       0       0       0       0       0',
        '       0       0       0       F       T       0       0       0       1       0       0       0',
        '       0       0       F       F       0       0       0       F       0       0       F',
        '       0',
        '       1       1       0       F       T       F       0       3      30       T       T       F       0',
        '       0     .01       F     .01      20      .3     .01       T       5       T       F',
        'Section - Expansion and Contraction Coefficients',
        '       0       0       2',
        'Section - Expansion and Contraction Coefficients (Unsteady)',
        '       2',
        '       0       0       0       0',
        'Section - Minor Loss Coefficients (Geometry)',
        '       2',
        '       0       0',
        'Section - Job Control',
        '       1       1       1       F       F       1',
        'Section - Junction Information',
        'Section - Reach Boundaries',
        '       T       T       2       3        Fake River             F       F',
        'Section - Encroachment Data',
        'Section - Observed Water Surface Data',
        'Section - Roughness Change Factors',
        'Section - Breach Data',
        'Section - Storage Area Data',
        'SA     8                                                        %s' % area_padded,
        '       0       0       0    %4d       T' % n_perim,
        'Section - Storage Area Connection Data',
        'Section - Pump Station Data',
        'Section - River Reach Data',
        'NODE   1Fake Reach      100                  100               0',
        '       1       F',
        '       4',
        '       0      10       0       0     100       0     100      10',
        'Section - XS Manning\'s/Roughness Data',
        '       1       F       F       0',
        '       0     .03',
        '       0       0       0       0       F       F       0     100       F       F       0       0       0       0',
        '       1       0       F',
        '       0     100       0',
        '       0       T',
        '       0       0       0',
        '      .5     .47      20',
        'NODE   1Fake Reach      0                                      0',
        '       1       F',
        '       4',
        '       0      10       0       0     100       0     100      10',
        'Section - XS Manning\'s/Roughness Data',
        '       1       F       F       0',
        '       0     .03',
        '       0       0       0       0       F       F       0     100       F       F       0       0       0       0',
        '       1       0       F',
        '       0     100       0',
        '       0       T',
        '       0       0       0',
        '      .5     .47      20',
    ]
    with open(x01_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    logger.info("  Generated .x01: %s (area=%s, perimeter points=%d)",
                x01_path.name, area_name, n_perim)


def _strip_results(src_hdf: Path, dst_hdf: Path) -> None:
    """Copy src_hdf → dst_hdf omitting Results and other result-only groups.

    This mirrors the remove_HDF5_Results_Sed.py utility and is used by
    Workflow A to prepare a fresh p01.tmp.hdf from an existing p01.hdf.
    """
    SKIP = {
        'Results',
        'Bed Time Series',
        'DSS Time Series',
        'Transport Time Series',
        'Unsteady Time Series',
    }
    with h5py.File(str(src_hdf), 'r') as src, h5py.File(str(dst_hdf), 'w') as dst:
        # Copy root-level attributes
        for attr_key, attr_val in src.attrs.items():
            try:
                dst.attrs[attr_key] = attr_val
            except Exception:
                pass
        # Copy all groups/datasets except results groups
        for key in src.keys():
            if key not in SKIP:
                src.copy(key, dst)


def _g01hdf_has_tables(g01_hdf_path: Path) -> bool:
    """Return True if the g01.hdf already contains computed hydraulic tables."""
    try:
        with h5py.File(str(g01_hdf_path), 'r') as hf:
            geo = hf.get('Geometry/2D Flow Areas', {})
            for area_name in geo.keys():
                if area_name == 'Attributes':
                    continue
                area_grp = geo[area_name]
                if isinstance(area_grp, h5py.Group) and \
                        'Cells Volume Elevation Info' in area_grp.keys():
                    return True
    except Exception:
        pass
    return False


def _find_terrain_layer(project_dir: Path, g01_hdf: Path = None):
    """The terrain RASMapper uses: the geometry's own Terrain Filename (a
    terrain HDF listing its rasters) when it has one, else a terrain HDF in
    Terrain/, else a terrain GeoTIFF."""
    if g01_hdf is not None and g01_hdf.exists():
        with h5py.File(str(g01_hdf), "r") as hf:
            name = hf["Geometry"].attrs.get("Terrain Filename", b"").decode("utf-8", "replace")
        if name:
            path = (project_dir / name.replace("\\", "/")).resolve()
            if path.exists():
                return path
    for hdf in sorted((project_dir / "Terrain").glob("*.hdf")):
        with h5py.File(str(hdf), "r") as hf:
            if "Terrain" in hf:
                return hdf
    return _find_terrain_tif(project_dir)


def _find_layer(project_dir: Path, g01_hdf: Path, attr: str, default: str) -> Path:
    """A land classification layer the geometry names (its HDF): at the
    geometry's path, else anywhere in the project under the same file name,
    else Land Classification/<default>."""
    if g01_hdf.exists():
        with h5py.File(str(g01_hdf), "r") as hf:
            name = hf["Geometry"].attrs.get(attr, b"").decode("utf-8", "replace").replace("\\", "/")
        if name:
            path = (project_dir / name).resolve()
            if path.exists():
                return path
            found = sorted(project_dir.rglob(Path(name).name))
            if found:
                return found[0]
    return project_dir / "Land Classification" / default


def _find_terrain_tif(project_dir):
    """Find terrain GeoTIFF in project directory."""
    terrain_dir = project_dir / "Terrain"
    if terrain_dir.exists():
        tifs = list(terrain_dir.glob("*.tif"))
        if tifs:
            return tifs[0]
    # Try project dir directly
    tifs = list(project_dir.glob("Terrain*.tif"))
    return tifs[0] if tifs else None


if __name__ == "__main__":
    main()
