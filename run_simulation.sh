#!/usr/bin/env bash
# =============================================================================
# run_simulation.sh — HEC-RAS 2D on Linux, zip in → zip out
#
# Usage:
#   ./run_simulation.sh <input.zip> [output.zip]
#
# Environment:
#   HECRAS_VERSION   engine version under bin/ (default 6.6; also 7.0.1)
#   HECRAS_RUNNER    native | docker (default: native on Linux x86_64)
#   EVENT_CONDITIONS path to an HDF whose boundary conditions Workflows B/C
#                    should use (see ras_preprocess.py --event-conditions)
#
# The input zip holds one HEC-RAS project folder:
#
#   MyProject.zip
#   └── MyProject/
#       ├── MyProject.g01          geometry text (seed points, perimeter, BC lines)
#       ├── MyProject.g01.hdf      mesh (Workflows A and B; absent for C)
#       ├── MyProject.p01          plan text
#       ├── MyProject.p01.hdf      GUI results (Workflow A only)
#       ├── MyProject.b01 / .x01   solver inputs (written if missing)
#       ├── MyProject.u01          unsteady flow file
#       ├── Terrain/*.tif          terrain (required)
#       └── Land Classification/   optional: LandCover.hdf + .tif, Infiltration.hdf
#
# The output zip holds <name>.g01.hdf, <name>.p01.hdf (results) and the
# computation log <name>.bco01.
#
# Requires Python 3 with requirements.txt installed, unzip and zip; plus
# Docker with the hecras-linux:<version> image when HECRAS_RUNNER=docker.
# =============================================================================
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <input.zip> [output.zip]" >&2
    exit 1
fi

_abspath() {
    case "$1" in
        /*) echo "$1" ;;
        *)  echo "$PWD/$1" ;;
    esac
}

INPUT_ZIP="$(_abspath "$1")"
OUTPUT_ZIP="$(_abspath "${2:-$(basename "${1%.zip}")_results.zip}")"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
HECRAS_VERSION="${HECRAS_VERSION:-6.6}"
PY="${PYTHON:-python3}"

echo "HEC-RAS ${HECRAS_VERSION} Linux runner"
echo "  Input  : $INPUT_ZIP"
echo "  Output : $OUTPUT_ZIP"

for cmd in "$PY" unzip zip; do
    command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: '$cmd' not found in PATH" >&2; exit 1; }
done

WORKDIR="$(mktemp -d)"
trap 'rm -rf "$WORKDIR"' EXIT

echo "[1/4] Extracting ..."
unzip -q "$INPUT_ZIP" -d "$WORKDIR/input"

# Project folder = the folder holding the first .g01 / .g01.hdf (skip __MACOSX)
G01="$(find "$WORKDIR/input" \( -name '*.g01' -o -name '*.g01.hdf' \) ! -path '*/__MACOSX/*' | sort | head -1)"
if [[ -z "$G01" ]]; then
    echo "ERROR: no .g01 or .g01.hdf file in $INPUT_ZIP" >&2
    exit 1
fi
PROJECT_DIR="$(dirname "$G01")"
echo "  Project folder: ${PROJECT_DIR#"$WORKDIR/input/"}"

RUN_DIR="$WORKDIR/run"
mkdir -p "$RUN_DIR"

echo "[2/4] Preprocessing (workflow detected automatically) ..."
PRE_ARGS=("$PROJECT_DIR" --output-dir "$RUN_DIR")
[[ -n "${EVENT_CONDITIONS:-}" ]] && PRE_ARGS+=(--event-conditions "$(_abspath "$EVENT_CONDITIONS")")
"$PY" "$SCRIPT_DIR/ras_preprocess.py" "${PRE_ARGS[@]}"

G01_HDF="$(ls "$RUN_DIR"/*.g01.hdf 2>/dev/null | head -1)"
PROJECT_NAME="$(basename "$G01_HDF" .g01.hdf)"

# Prefer the project's own solver inputs over generated ones
for f in "$PROJECT_DIR/$PROJECT_NAME".b01 "$PROJECT_DIR/$PROJECT_NAME".x01; do
    [[ -e "$f" ]] && cp "$f" "$RUN_DIR/"
done

echo "[3/4] Running RasGeomPreprocess + RasUnsteady ..."
"$PY" "$SCRIPT_DIR/hecras_run.py" "$RUN_DIR" "$PROJECT_NAME" --version "$HECRAS_VERSION" \
    ${HECRAS_RUNNER:+--runner "$HECRAS_RUNNER"}

echo "[4/4] Packaging ..."
OUTPUT_FILES=()
for pat in "*.g01.hdf" "*.p01.hdf" "*.bco01"; do
    for f in "$RUN_DIR"/$pat; do
        [[ -e "$f" ]] && OUTPUT_FILES+=("$(basename "$f")")
    done
done
rm -f "$OUTPUT_ZIP"
(cd "$RUN_DIR" && zip -q "$OUTPUT_ZIP" "${OUTPUT_FILES[@]}")
echo "Done: $OUTPUT_ZIP (${OUTPUT_FILES[*]})"
