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
MULTI_N = "Chippewa_2D (Multiple Face Mann n)"
PRODUCTION = "production model (calibration region, terrain modification)"

PORTED, PARTIAL, MISSING, SUBST, PRIM, NA = (
    "ported", "partly ported", "not ported", "substituted", "primitive", "not applicable")

# (class regex, method regex, status, hecras_linux counterpart, verified by / note)
RULES = [
    # -- the driver -------------------------------------------------------------
    (r"RASD2FlowArea$", r"^CreatePropertyTables$", PARTIAL, "preprocess.compute_face_properties, compute_cell_properties",
     "2D bridge profiles not ported"),
    (r"RASD2FlowArea/.*", r"^_Lambda\$__1[34]$", PORTED, "preprocess._composite_cells", "Chippewa_2D (Manning's n)"),
    (r"RASD2FlowArea$", r"^ComputeAssociationCellAverage", PORTED, "preprocess._composite_cells", "Chippewa_2D (Manning's n)"),
    (r"RASD2FlowArea$", r"^ProcessLandCoverPropertyTables$", PORTED, "preprocess.sample_land_cover", EXAMPLES + ", Chippewa_2D"),
    (r"RASD2FlowArea$", r"^ProcessClassificationPropertyTables$", PORTED, "preprocess.sample_infiltration",
     "BEC, Muncie, AFTER_RUN (cell averages unverified)"),
    (r"RASD2FlowArea$", r"^(Effective|Get)(CellAverageLC|MultipleFaceMannN)$", PORTED, "preprocess.table_settings", "Chippewa_2D"),
    (r"RASD2FlowArea$", r"^Get(FilterTol|FilterCell|FaceLaminar|FaceAreaConveyance|ManningsN)", PORTED, "preprocess.table_settings", EXAMPLES),
    (r"RASD2FlowArea", r"(Dataset|get_Geometry|\.c?ctor)$", NA, "", "names and object setup"),
    (r"RASD2FlowArea/", r"^_Lambda", NA, "", "delegate glue (progress, writers)"),
    # -- terrain ----------------------------------------------------------------
    (r"TerrainLayer(/.*)?$", r"^(ComputeProfilesOnRaster|ComputeProfilesOnTile|ComputePointElevationOnTile)$", PORTED,
     "property_tables._raster_profile, _tile_profile", EXAMPLES),
    (r"TerrainLayer$", r"^ComputeProfiles$", PARTIAL, "property_tables.face_profile", "faces crossing terrain modifications (ObstructProfiles) refused"),
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
    # The stitch TIN: CreatePropertyTables builds it from stitch segments only, so
    # perimeter, breakline, mass-point and quad-injection branches never run.
    (r"RasMapperLib\.DelaunayTINFactory(/.*)?$", r"^(InsertPerimeterPoints|InsertBreaklinePoints|InsertMassPoints|DelaunayInsertPoints|EnforceBreakLineEdges|EnforcePerimeterEdges|EnforceEdgeSegments|EnforceInternalMaxEdgeLengthConstraint|MarkBadTrianglesAsDeleted_MultiPart|CheckIfCurTriangleAndNeighborFormQuad|DeleteTrianglesLeftOfStitchEdge_DeleteAdjacentTriangles|TriangleCenter|TriangleMaxEdgeLength|TriangleArea|TrianglesMaxM|AddSuperTriangle|IsColinear|TrianglePoint(Before|After)Point|WriteInputDataToImmediateWindow|Case123|ResetGlobalVars)$",
     NA, "", "TIN branches unused for stitches (perimeter, breaklines, mass points)"),
    (r"RasMapperLib\.TINFactory$", r"^(Validate|FindTriangleContainingPoint2|e012After|Debug.*|Export.*)$", NA, "", "TIN branches unused for stitches"),
    (r"RasMapperLib\.DevelopmentTin(/.*)?$", r"^(Rasterize|RemoveDeletedTriangles|\.ctor|get_Triangles|TrianglePointMs)$", PORTED, "tin.StitchTin.compute, tin.paint", "BaldEagle g09, VA, Chippewa_2D"),
    (r"RasMapperLib\.DevelopmentTin(/.*)?$", r".", NA, "", "point bookkeeping after removal (no effect on painted values)"),
    (r"RasMapperLib\.(DelaunayTINFactory|TINFactory|TinHelpers\..*|PlaneM|CircleM)(/.*)?$", r".", PORTED,
     "tin.StitchTin", "BaldEagle g09, VA, Chippewa_2D"),
    (r"RasMapperLib\.RasterizeTriangle$", r"^ComputeCellsZ$", PORTED, "tin.paint", "BaldEagle g09, VA, Chippewa_2D"),
    (r"RasMapperLib\.TriangleShape$", r".", PORTED, "property_tables.TerrainLayer.stitch_profile", "BaldEagle g09"),
    # -- faces and cells ------------------------------------------------------------
    (r"ProfileHydraulicPropertiesTable(/.*)?$", r"^\.ctor$", None, "", ""),  # decided by size below
    (r"ProfileHydraulicPropertiesTable(/.*)?$", r"^ToString$", NA, "", "debug text"),
    (r"ProfileHydraulicPropertiesTable(/.*)?$", r".", PORTED, "property_tables.face_table, _HP, centroid_at_area", EXAMPLES),
    (r"RasMapperLib\.Profile$", r"^(GetSubSection|GetSubsectionInternal)$", PORTED, "property_tables._sub_section", MULTI_N),
    (r"RasMapperLib\.Profile$", r"^Average$", MISSING, "", "face variable averages (Composite LC)"),
    (r"RasMapperLib\.Profile$", r".", PORTED, "property_tables.face_profile", EXAMPLES),
    (r"RasMapperLib\.Point2DCollection$", r"^ToString$", NA, "", "debug text"),
    (r"RasMapperLib\.Point2DCollection$", r".", PORTED, "property_tables._dp_keep, _dotnet_sort", EXAMPLES),
    (r"RasMapperLib\.Histogram$", r".", PORTED, "property_tables.cell_table", CELLS_ALSO),
    (r"RasMapperLib\.MeshFV2D(/.*)?$", r"^(CellPolygonPoints|CellPolygon|ComputeFaceProfiles|FaceSegment|FacePolyline|FacePolylineLength|PointsOnFace|CellIsVirtual|get_NonVirtualCellCount|Set(Cell|Face)MinElevation|Cell|get_(Cell|Face|FacePoint)Count)$",
     PORTED, "preprocess._cell_polygon, _face_polyline, read_mesh_from_g01hdf", EXAMPLES),
    (r"RasMapperLib\.MeshFV2D(/.*)?$", r".", NA, "", "mesh bookkeeping (extents, messages, names)"),
    (r"RasMapperLib\.RasterM$", r".", PORTED, "property_tables._cell_raster, Terrain.locate", EXAMPLES),
    # -- land cover, Manning's n, infiltration ------------------------------------------
    # Closure numbers are RASGeometry's methods: 665 ExtractLCVariablePolylines, 666 ExtractManningsNPolylines,
    # 667 raster extraction, 668 vector extraction, 669 combine, 670 fill, 671 sort, 672-673 tiles.
    (r"RASGeometry/_Closure\$__$", r"^\.c?ctor$", NA, "", "object setup"),
    (r"RASGeometry/_Closure\$__(665-\d)?$", r"^(_Lambda\$__665-|_Lambda|\.ctor)", None, "", ""),
    (r"RASGeometry/_Closure\$__(66[89]-\d)?$", r"^(_Lambda\$__66[89]-|_Lambda|\.ctor)", None, "", ""),
    (r"RASGeometry/_Closure\$__(67[0-3]-\d)?$", r".", PORTED, "property_tables.face_n_profiles", MULTI_N),
    (r"RASGeometry$", r"^(ExtractManningsNPolylines|ExtractLCVariablePolylinesRegionsRasterExtraction|ExtractLCVariablePolylineFillNoData|ExtractLCVariablePolyineSort|ExtractLCVariableSampleTile|ExtractLCVariableTilesToProcessFromTileRaster)$",
     PORTED, "property_tables.face_n_profiles", MULTI_N),
    (r"RASGeometry$", r"^ExtractLCVariableSubsampleTile$", NA, "", "only with sub-sampling, which no table caller turns on"),
    (r"RASGeometry$", r"^(ExtractLCVariablePolylinesRegionsVectorExtraction|ExtractLCVariablePolylineRegionsCombine)$",
     PORTED, "property_tables._n_from_regions, _combine_n, _SegmentIndex", PRODUCTION),
    (r"RASGeometry$", r"^ExtractLCVariablePolylines$", MISSING, "", "face variable averages (Composite LC)"),
    (r"RasMapperLib\.FinalNValueLayer$", r".", PORTED, "preprocess.sample_land_cover", EXAMPLES),
    (r"RasMapperLib\.RasLandCoverRegions(/.*)?$", r"^MD5", NA, "", "change detection"),
    (r"RasMapperLib\.RasLandCoverRegions(/.*)?$", r".", PORTED, "preprocess.read_calibration", PRODUCTION),
    (r"RasMapperLib\.LandCoverLayer(/.*)?$", r"^(OverridePolygonFiltered|GetPolysAsMappedTuples|get_OverridePolygonLayer)$", MISSING, "", "override regions"),
    (r"RasMapperLib\.LandCoverLayer(/.*)?$", r"^(get_Resampler|Rasterize|GetTileValues|get_Classification|GetVirtualClassification|EvaluateVirtualDependencies|get_VirtualDependencyMode|get_Parameters|get_LCType|GetClassificationName)$",
     PORTED, "preprocess._classify_points, sample_infiltration", "BEC, Muncie, AFTER_RUN"),
    (r"RasMapperLib\.LandCoverLayer(/.*)?$", r".", NA, "", "layer plumbing (names, validation, files)"),
    (r"RasMapperLib\.LandCoverLayerHelper$", r".", PORTED, "preprocess.sample_infiltration (Variables table)", "BEC, Muncie, AFTER_RUN"),
    (r"RasMapperLib\.LandCoverCategor(ies|y)$", r".", NA, "", "default class lists for new layers"),
    (r"Classifications\.ParameterResampler(/.*)?$", r"^(Resample|<Resample>)", PORTED, "preprocess._composite_cells", "Chippewa_2D (Manning's n)"),
    (r"Classifications\.ParameterResampler(/.*)?$", r"(PolygonVariableMap|Polygons(Highest|Lowest))", PORTED,
     "preprocess.sample_land_cover (regions at points)", PRODUCTION),
    (r"Classifications\.ParameterResampler(/.*)?$", r".", PORTED, "preprocess.sample_land_cover, sample_infiltration", "BEC, Muncie, AFTER_RUN"),
    (r"Classifications\.(ClassificationResampler|MergedClassificationResampler)(/.*)?$", r"^(Resample|<Resample>)", PORTED,
     "preprocess._ClassRaster.classes_at", "Chippewa_2D (Manning's n)"),
    (r"Classifications\.", r".", PORTED, "preprocess._classify_points (virtual layers)", "BEC, Muncie, AFTER_RUN"),
    (r"^Geospatial\.Vectors\.Rasterizers\.", r".", MISSING, "", "override regions (rasterized onto a cell's raster)"),
    (r"^Geospatial\.Vectors\.(DelaunayTinFactory|Tin|TinFactory|Triangle.*|EdgeDefinition|TriangulationParameters|Circle)(/.*)?$", r".", PORTED,
     "tin.PolygonTin", PRODUCTION),
    (r"^Geospatial\.Rasters\.", r".", PORTED, "preprocess._ClassRaster (nearest pixel, pyramid level)", "Chippewa_2D (Manning's n)"),
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
            if status is None:
                if crx.startswith("RASGeometry"):   # closures, by the method they belong to
                    num = re.search(r"(665|668|669)-", cls + method)
                    if num and num.group(1) == "665":
                        return MISSING, "", "face variable averages (Composite LC)"
                    if num:
                        return PORTED, "property_tables._n_from_regions, _combine_n", PRODUCTION
                    continue
                if il > 2000:   # ProfileHydraulicPropertiesTable constructors
                    return PORTED, "property_tables.face_table_n_profile", MULTI_N
                return PORTED, "property_tables.face_table", EXAMPLES
            return status, port, note
    return "unclassified", "", ""


PRIM_USED, PRIM_UNUSED = "primitive, used by ported code", "primitive, only for unported features"


def main(reach, stopped, md_out, tsv_out):
    rows = list(csv.DictReader(open(reach), delimiter="\t"))
    out = []
    key_of = {}
    for r in rows:
        cls = r["type"]
        method = re.sub(r"\(.*", "", r["method"])
        il = int(r["il_bytes"])
        status, port, note = classify(cls, method, il)
        key_of[f'{cls}::{r["method"]}'] = len(out)
        out.append([r["assembly"], cls, method, il, status, port, note])
    # A primitive is covered when ported code reaches it, directly or through
    # other primitives and not-applicable helpers.
    callees = defaultdict(list)
    for line in open(reach + ".edges"):
        a, b = line.rstrip("\n").split("\t")
        if a in key_of and b in key_of:
            callees[key_of[a]].append(key_of[b])
    stack = [i for i, o in enumerate(out) if o[4] in (PORTED, PARTIAL)]
    reached = set(stack)
    while stack:
        for j in callees[stack.pop()]:
            if j not in reached and out[j][4] in (PRIM, NA):
                reached.add(j)
                stack.append(j)
    for i, o in enumerate(out):
        if o[4] == PRIM:
            o[4] = PRIM_USED if i in reached else PRIM_UNUSED
            o[6] = "" if i in reached else "only reached from unported features"

    with open(tsv_out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["assembly", "type", "method", "il_bytes", "status", "hecras_linux", "verified_by_or_note"])
        for row in sorted(out, key=lambda x: (x[4], x[1], x[2])):
            w.writerow(row)

    tot = defaultdict(lambda: [0, 0])
    for _, _, _, il, status, _, _ in out:
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
        PRIM_USED: "geometry and math building blocks the ported code uses; reproduced inline",
        PRIM_UNUSED: "building blocks reached only from unported features",
        NA: "can't change table values (caches, names, debug text, setup)",
        "unclassified": "no rule yet",
    }
    for s in [PORTED, PARTIAL, MISSING, SUBST, PRIM_USED, PRIM_UNUSED, NA, "unclassified"]:
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
    labels = {"terrainmod": "Terrain modifications (polygons with one elevation are applied to cell "
                            "rasters, in property_tables._apply_modifications; other types, and faces "
                            "crossing any modification, are refused)",
              "feature": "Other features (not ported)", "model": "Project and layer loading (read from the HDF instead)",
              "io": "File I/O (h5py and GDAL instead)", "ui": "User interface", "progress": "Progress reporting"}
    # A call graph from before terrain modifications had their own category files them as features.
    obstruct = [n for n in st.get("feature", []) if "Obstruct" in n]
    st["feature"] = [n for n in st.get("feature", []) if "Obstruct" not in n]
    st["terrainmod"] = st.get("terrainmod", []) + obstruct
    for cat in ["terrainmod", "feature", "model", "io", "ui", "progress"]:
        names = sorted(set(st.get(cat, [])))
        L.append(f"- **{labels[cat]}:** {len(names)} entry points" +
                 (": " + ", ".join(f"`{n}`" for n in names) if cat in ("feature", "terrainmod") else ""))
    L.append("")
    L.append(f"Every method, with its status and counterpart: [`rasmapper-coverage.tsv`]({tsv_out.split('/')[-1]}).\n")
    open(md_out, "w").write("\n".join(L))
    print(f"{len(out)} methods; " + ", ".join(f"{s}: {tot[s][0]}" for s in tot))


if __name__ == "__main__":
    main(*sys.argv[1:5])
