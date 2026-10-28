#!/usr/bin/env python3
"""How much of RASMapper's property-table code hecras_linux covers.

    python3 tools/rasmapper_coverage.py reach.tsv reach.tsv.stopped \\
        docs/rasmapper-coverage.md docs/rasmapper-coverage.tsv

reach.tsv comes from tools/rasmapper_callgraph (every method reachable from
RASD2FlowArea.CreatePropertyTables in HEC's DLLs). Each method is classified
by the first matching rule below; a method no rule matches is reported as
unclassified, so the inventory can't silently miss one.
"""
from __future__ import annotations

import csv
import re
import sys
from collections import defaultdict

# Models whose GUI-computed tables (RASMapper 7.0.2 output) the port matches
# bit for bit; see README "Workflow B accuracy".
EXAMPLES = "5 examples, Weise_2D, BaldEagle g09"
CELLS_ALSO = EXAMPLES + ", Chippewa_2D (cells)"

PORTED, PARTIAL, MISSING, SUBST, PRIM, NA = (
    "ported", "partly ported", "not ported", "substituted", "primitive", "not applicable")

# (class regex, method regex, status, hecras_linux counterpart, verified by / note)
RULES = [
    # -- the driver -------------------------------------------------------------
    (r"RASD2FlowArea$", r"^CreatePropertyTables$", PARTIAL, "preprocess.compute_face_properties, compute_cell_properties",
     "Composite LC and Multiple Face Mann n branches refused"),
    (r"RASD2FlowArea/.*", r"^_Lambda\$__1[34]$", MISSING, "", "Composite LC: per-cell association averages (parallel loop)"),
    (r"RASD2FlowArea$", r"^ComputeAssociationCellAverage", MISSING, "", "Composite LC"),
    (r"RASD2FlowArea$", r"^ProcessLandCoverPropertyTables$", PARTIAL, "preprocess.sample_land_cover",
     "point sampling ported; cell averaging (Composite LC) not"),
    (r"RASD2FlowArea$", r"^ProcessClassificationPropertyTables$", PARTIAL, "preprocess.sample_infiltration, _classify_points",
     "point sampling ported; override regions and cell averaging not"),
    (r"RASD2FlowArea$", r"^(Effective|Get)(CellAverageLC|MultipleFaceMannN)$", PORTED, "preprocess.table_settings", "refused with an error"),
    (r"RASD2FlowArea$", r"^Get(FilterTol|FilterCell|FaceLaminar|FaceAreaConveyance|ManningsN)", PORTED, "preprocess.table_settings", EXAMPLES),
    (r"RASD2FlowArea", r"(Dataset|get_Geometry|\.c?ctor)$", NA, "", "names and object setup"),
    (r"RASD2FlowArea/", r"^_Lambda", NA, "", "delegate glue (progress, writers)"),
    # -- terrain ----------------------------------------------------------------
    (r"TerrainLayer(/.*)?$", r"^(ComputeProfilesOnRaster|ComputeProfilesOnTile|ComputePointElevationOnTile)$", PORTED,
     "property_tables._raster_profile, _tile_profile", EXAMPLES),
    (r"TerrainLayer$", r"^ComputeProfiles$", PARTIAL, "property_tables.face_profile", "terrain modifications (ObstructProfiles) not ported"),
    (r"TerrainLayer$", r"^(ResampleInternal|Resample|BufferRenderBlockData)$", PORTED, "property_tables.cell_table, Terrain.sample", CELLS_ALSO),
    (r"TerrainLayer$", r"^Level0StitchTIN$", PORTED, "property_tables.TerrainLayer.stitch_profile", "BaldEagle g09"),
    (r"TerrainLayer$", r"^ClearLRUTilesConcurrent$", NA, "", "tile cache"),
    (r"TerrainLayer$", r"^get_", PORTED, "property_tables.TerrainLayer", EXAMPLES),
    (r"TerrainLayer/", r".", PORTED, "property_tables (profile and resample loops)", EXAMPLES),
    (r"TerrainLayer$", r"\.c?ctor$", NA, "", "object setup"),
    (r"RasMapperLib\.Terrain\.", r".", PORTED, "property_tables.Terrain, _interpolate", EXAMPLES),
    (r"RasMapperLib\.(RasterizePolygon|RasterizeSegment|RasterizeTriangle)(/.*)?$", r".", PORTED,
     "property_tables._rasterize, _walk_cells, _paint_triangles", EXAMPLES),
    (r"RasMapperLib\.Tin(/.*)?$", r"^(ComputeProfiles|_Lambda\$__[12]|TriangleShape|TriangleEdgeSegment|TriangleCount|TriangleArea)$", PORTED,
     "property_tables.TerrainLayer.stitch_profile", "BaldEagle g09"),
    (r"RasMapperLib\.Tin(/.*)?$", r".", PORTED, "property_tables.TerrainLayer (bounding-box search)", "BaldEagle g09"),
    (r"RasMapperLib\.(DelaunayTINFactory|TINFactory|DevelopmentTin|TinHelpers\..*|PlaneM|CircleM)(/.*)?$", r".", SUBST,
     "property_tables._stitch_triangles (scipy Delaunay)", "stitch TIN for cells; RASMapper's own TIN to be ported"),
    (r"RasMapperLib\.TriangleShape$", r".", PORTED, "property_tables.TerrainLayer.stitch_profile", "BaldEagle g09"),
    # -- faces and cells ------------------------------------------------------------
    (r"ProfileHydraulicPropertiesTable(/.*)?$", r"^\.ctor$", None, "", ""),  # decided by size below
    (r"ProfileHydraulicPropertiesTable(/.*)?$", r"^ToString$", NA, "", "debug text"),
    (r"ProfileHydraulicPropertiesTable(/.*)?$", r".", PORTED, "property_tables.face_table, _HP, centroid_at_area", EXAMPLES),
    (r"RasMapperLib\.Profile$", r"^(GetSubSection|GetSubsectionInternal|Average)$", MISSING, "", "Multiple Face Mann n subsections"),
    (r"RasMapperLib\.Profile$", r".", PORTED, "property_tables.face_profile", EXAMPLES),
    (r"RasMapperLib\.Point2DCollection$", r"^ToString$", NA, "", "debug text"),
    (r"RasMapperLib\.Point2DCollection$", r".", PORTED, "property_tables._dp_keep, _dotnet_sort", EXAMPLES),
    (r"RasMapperLib\.Histogram$", r".", PORTED, "property_tables.cell_table", CELLS_ALSO),
    (r"RasMapperLib\.MeshFV2D(/.*)?$", r"^(CellPolygonPoints|CellPolygon|ComputeFaceProfiles|FaceSegment|FacePolyline|FacePolylineLength|PointsOnFace|CellIsVirtual|get_NonVirtualCellCount|Set(Cell|Face)MinElevation|Cell|get_(Cell|Face|FacePoint)Count)$",
     PORTED, "preprocess._cell_polygon, _face_polyline, read_mesh_from_g01hdf", EXAMPLES),
    (r"RasMapperLib\.MeshFV2D(/.*)?$", r".", NA, "", "mesh bookkeeping (extents, messages, names)"),
    (r"RasMapperLib\.RasterM$", r".", PORTED, "property_tables._cell_raster, Terrain.locate", EXAMPLES),
    # -- land cover, Manning's n, infiltration ------------------------------------------
    (r"RASGeometry(/.*)?$", r"^(Extract|_Lambda|\.c?ctor)", MISSING, "", "Multiple Face Mann n: n profiles along faces"),
    (r"RasMapperLib\.FinalNValueLayer$", r".", PORTED, "preprocess.sample_land_cover", EXAMPLES),
    (r"RasMapperLib\.RasLandCoverRegions(/.*)?$", r"^MD5", NA, "", "change detection"),
    (r"RasMapperLib\.RasLandCoverRegions(/.*)?$", r".", MISSING, "", "calibration regions (Manning's n polygons)"),
    (r"RasMapperLib\.LandCoverLayer(/.*)?$", r"^(OverridePolygonFiltered|GetPolysAsMappedTuples|get_OverridePolygonLayer)$", MISSING, "", "override regions"),
    (r"RasMapperLib\.LandCoverLayer(/.*)?$", r"^(get_Resampler|Rasterize|GetTileValues|get_Classification|GetVirtualClassification|EvaluateVirtualDependencies|get_VirtualDependencyMode|get_Parameters|get_LCType|GetClassificationName)$",
     PORTED, "preprocess._classify_points, sample_infiltration", "BEC, Muncie, AFTER_RUN"),
    (r"RasMapperLib\.LandCoverLayer(/.*)?$", r".", NA, "", "layer plumbing (names, validation, files)"),
    (r"RasMapperLib\.LandCoverLayerHelper$", r".", PORTED, "preprocess.sample_infiltration (Variables table)", "BEC, Muncie, AFTER_RUN"),
    (r"RasMapperLib\.LandCoverCategor(ies|y)$", r".", NA, "", "default class lists for new layers"),
    (r"Classifications\.ParameterResampler(/.*)?$", r"^(Resample|<Resample>)", MISSING, "", "Composite LC: resampling onto a cell's raster"),
    (r"Classifications\.ParameterResampler(/.*)?$", r"(PolygonVariableMap|Polygons(Highest|Lowest))", MISSING, "", "override regions"),
    (r"Classifications\.ParameterResampler(/.*)?$", r".", PORTED, "preprocess.sample_land_cover, sample_infiltration", "BEC, Muncie, AFTER_RUN"),
    (r"Classifications\.(ClassificationResampler|MergedClassificationResampler)(/.*)?$", r"^(Resample|<Resample>)", MISSING, "", "Composite LC"),
    (r"Classifications\.", r".", PORTED, "preprocess._classify_points (virtual layers)", "BEC, Muncie, AFTER_RUN"),
    (r"^Geospatial\.Vectors\.Rasterizers\.", r".", MISSING, "", "Composite LC / override regions"),
    (r"^Geospatial\.Vectors\.(DelaunayTinFactory|Tin|TinFactory|Triangle.*|EdgeDefinition|TriangulationParameters|Circle)(/.*)?$", r".", MISSING, "", "override regions (polygon containment TINs)"),
    (r"^Geospatial\.Rasters\.", r".", MISSING, "", "Composite LC (raster resampling of classification layers)"),
    (r"^Geospatial\.IO\.", r".", NA, "", "TIFF reading (GDAL instead)"),
    # -- shared building blocks ----------------------------------------------------------
    (r"^RasMapperLib\.(SegmentM|Extent|PointM|PointMs|Polygon|Polyline|MultiPoint|Point2D|Feature|SpatialIndex(`1)?|DistinctSpatialIndex|PointOnPolyline|PointIntersect|Utilities\.(Vector|Math))(/.*)?$", r".", PRIM,
     "numpy / inline arithmetic", ""),
    (r"^Geospatial\.(Vectors|DataStructures)\.", r".", PRIM, "numpy / inline arithmetic", ""),
    (r"^Utility\.Math$", r".", PRIM, "numpy", ""),
    (r"^(Utility\.|RasMapperLib\.(SharedData|Utilities\.(Path|DSS|Hash|Converter|General)|DataTableExtensions|DictionaryExtensionMethods|ListExtensionMethods|ImageLayer|RASMapper|MeshFeatureLayer|ManningsNPolygonLayer|HDFBufferedWriter`1)\b)", r".", NA, "",
     "infrastructure (caches, threading, paths, HDF writing)"),
]


def classify(cls: str, method: str, il: int):
    for crx, mrx, status, port, note in RULES:
        if re.search(crx, cls) and re.search(mrx, method):
            if status is None:   # ProfileHydraulicPropertiesTable constructors
                if il > 2000:
                    return MISSING, "", "Multiple Face Mann n constructor (subsections)"
                return PORTED, "property_tables.face_table", EXAMPLES
            return status, port, note
    return "unclassified", "", ""


def main(reach, stopped, md_out, tsv_out):
    rows = list(csv.DictReader(open(reach), delimiter="\t"))
    out = []
    for r in rows:
        cls = r["type"]
        method = re.sub(r"\(.*", "", r["method"])
        il = int(r["il_bytes"])
        status, port, note = classify(cls, method, il)
        out.append((r["assembly"], cls, method, il, status, port, note))

    with open(tsv_out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["assembly", "type", "method", "il_bytes", "status", "hecras_linux", "verified_by_or_note"])
        for row in sorted(out, key=lambda x: (x[4], x[1], x[2])):
            w.writerow(row)

    tot = defaultdict(lambda: [0, 0])
    for *_, il, status, _, _ in [(a, c, m, i, s, p, n) for a, c, m, i, s, p, n in out]:
        tot[status][0] += 1
        tot[status][1] += il
    compute = [PORTED, PARTIAL, MISSING, SUBST]
    c_il = sum(tot[s][1] for s in compute)
    gaps = defaultdict(lambda: [0, 0])
    for a, c, m, il, s, p, n in out:
        if s in (MISSING, SUBST, PARTIAL):
            gaps[n][0] += 1
            gaps[n][1] += il

    st_rows = list(csv.DictReader(open(stopped), delimiter="\t"))
    st = defaultdict(list)
    for r in st_rows:
        st[r["category"]].append(f'{r["type"]}::{r["method"]}')

    L = []
    L.append("# RASMapper property-table coverage\n")
    L.append("How much of RASMapper's own code for *Compute Property Tables* "
             "(HEC-RAS 7.0.2, `RASD2FlowArea.CreatePropertyTables`) `hecras_linux` "
             "reproduces. Generated by `tools/rasmapper_coverage.py` from the call graph "
             "that `tools/rasmapper_callgraph` reads from the IL of HEC's DLLs; regenerate "
             "both after porting more.\n")
    L.append("## Summary\n")
    L.append(f"{len(out)} methods are reachable from the table computation once user interface, "
             "progress reporting, file I/O, project loading and other features are cut off "
             "(those are listed at the end).\n")
    L.append("| Status | Methods | IL bytes | Meaning |")
    L.append("|---|---|---|---|")
    meaning = {
        PORTED: "reproduced; tables match RASMapper's bit for bit on the models named",
        PARTIAL: "the paths the models use are reproduced; named branches are not",
        MISSING: "not reproduced; the feature it serves is refused or absent",
        SUBST: "reproduced with a different algorithm (results matched so far, not guaranteed)",
        PRIM: "geometry and math building blocks, carried along with the methods above",
        NA: "can't change table values (caches, names, debug text, setup)",
        "unclassified": "no rule yet",
    }
    for s in [PORTED, PARTIAL, MISSING, SUBST, PRIM, NA, "unclassified"]:
        if tot[s][0]:
            L.append(f"| {s} | {tot[s][0]} | {tot[s][1]:,} | {meaning[s]} |")
    L.append("")
    full = tot[PORTED][1]
    L.append(f"Of the computational code (ported, partly ported, not ported and substituted: "
             f"{c_il:,} IL bytes), **{100 * full / c_il:.0f} %** is fully ported, "
             f"{100 * tot[PARTIAL][1] / c_il:.0f} % partly, {100 * tot[SUBST][1] / c_il:.0f} % substituted "
             f"and {100 * tot[MISSING][1] / c_il:.0f} % not ported.\n")
    L.append("## What is missing, by feature\n")
    L.append("| Feature | Methods | IL bytes |")
    L.append("|---|---|---|")
    for note, (n, il) in sorted(gaps.items(), key=lambda kv: -kv[1][1]):
        L.append(f"| {note} | {n} | {il:,} |")
    L.append("")
    L.append("## Cut off from the call graph\n")
    L.append("Reached from the table computation but not followed (`tools/rasmapper_callgraph/stop.tsv`):\n")
    labels = {"feature": "Other features (not ported)", "model": "Project and layer loading (read from the HDF instead)",
              "io": "File I/O (h5py and GDAL instead)", "ui": "User interface", "progress": "Progress reporting"}
    for cat in ["feature", "model", "io", "ui", "progress"]:
        names = sorted(set(st.get(cat, [])))
        L.append(f"- **{labels[cat]}:** {len(names)} entry points" +
                 (": " + ", ".join(f"`{n}`" for n in names) if cat == "feature" else ""))
    L.append("")
    L.append(f"Every method, with its status and counterpart: [`rasmapper-coverage.tsv`]({tsv_out.split('/')[-1]}).\n")
    open(md_out, "w").write("\n".join(L))
    print(f"{len(out)} methods; " + ", ".join(f"{s}: {tot[s][0]}" for s in tot))


if __name__ == "__main__":
    main(*sys.argv[1:5])
