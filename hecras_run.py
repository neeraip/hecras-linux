#!/usr/bin/env python3
"""Run the HEC-RAS Linux engines on a preprocessed plan.

    python3 hecras_run.py <work_dir> <project_name> [--plan 01]
                          [--version 6.6|7.0.1] [--runner native|docker]

<work_dir> must hold <project>.p01.tmp.hdf, <project>.g01.hdf and the .b01 /
.x01 files that ras_preprocess.py writes. RasGeomPreprocess and RasUnsteady
are run in turn and the plan is renamed to <project>.p01.hdf.

Runners
  native  runs bin/<version>/ directly with libs/ on LD_LIBRARY_PATH.
          Needs Linux on x86_64 (glibc 2.28 or newer).
  docker  runs the hecras-linux:<version> image (see Dockerfile). Use this on
          macOS, Windows or ARM.
The default is native where it can work, docker otherwise.

The engines exit with status 0 even when they stop on an input error, so a
run counts as successful only if RasUnsteady reports that it finished and the
plan HDF holds results. Otherwise EngineError is raised with the engine log.
"""
from __future__ import annotations

import argparse
import os
import platform
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

REPO = Path(__file__).resolve().parent
BIN_ROOT = REPO / "bin"
LIBS = REPO / "libs"
DEFAULT_VERSION = "6.6"
DOCKER_IMAGE = "hecras-linux:{version}"

_VOLUME_ERR = re.compile(r"Overall Volume Accounting Error as percentage:\s*([-\d.eE+]+)")
_FINISHED = "Finished Unsteady Flow Simulation"


class EngineError(RuntimeError):
    """An engine stopped without producing results."""


@dataclass
class RunResult:
    plan_hdf: Path
    log: str
    volume_error_pct: float | None


def available_versions() -> list[str]:
    """Engine versions whose binaries are present under bin/."""
    if not BIN_ROOT.is_dir():
        return []
    return sorted(d.name for d in BIN_ROOT.iterdir()
                  if (d / "RasUnsteady").is_file() and (d / "RasUnsteady").stat().st_size > 1024)


def native_supported() -> bool:
    return platform.system() == "Linux" and platform.machine() in ("x86_64", "AMD64")


def docker_image_present(version: str) -> bool:
    if not shutil.which("docker"):
        return False
    r = subprocess.run(["docker", "image", "inspect", DOCKER_IMAGE.format(version=version)],
                       capture_output=True)
    return r.returncode == 0


def default_runner() -> str:
    return "native" if native_supported() else "docker"


def _engine_error(text: str) -> str | None:
    """First error block the engine printed, if any."""
    lines = [l.strip() for l in text.splitlines()]
    for i, line in enumerate(lines):
        if line.startswith("FONT= RED") or "HDF_ERROR" in line:
            block = [l for l in lines[i:i + 12] if l and not l.startswith("FONT=")]
            return " | ".join(block)
    return None


def run_engine(work_dir: Path, project: str, plan: str = "01",
               version: str = DEFAULT_VERSION, runner: str | None = None,
               timeout: int = 6 * 3600) -> RunResult:
    """Run RasGeomPreprocess + RasUnsteady; return the result plan and log."""
    work_dir = Path(work_dir).resolve()
    runner = runner or default_runner()
    plan = plan.zfill(2)
    tmp_hdf = f"{project}.p{plan}.tmp.hdf"
    out_hdf = work_dir / f"{project}.p{plan}.hdf"
    if not (work_dir / tmp_hdf).is_file():
        raise FileNotFoundError(work_dir / tmp_hdf)
    steps = [["RasGeomPreprocess", tmp_hdf, f"x{plan}"],
             ["RasUnsteady", tmp_hdf, f"x{plan}"]]

    log = ""
    if runner == "native":
        if not native_supported():
            raise EngineError("The native runner needs Linux on x86_64; use --runner docker")
        bin_dir = BIN_ROOT / version
        if not (bin_dir / "RasUnsteady").is_file():
            raise EngineError(f"No HEC-RAS {version} engines in {bin_dir} "
                              f"(installed: {', '.join(available_versions()) or 'none'})")
        env = dict(os.environ)
        env["LD_LIBRARY_PATH"] = os.pathsep.join(
            [str(LIBS), str(LIBS / "rhel_8"), str(LIBS / "mkl"), env.get("LD_LIBRARY_PATH", "")])
        for exe, *args in steps:
            try:
                p = subprocess.run([str(bin_dir / exe), *args], cwd=work_dir, env=env,
                                   capture_output=True, text=True, timeout=timeout)
            except subprocess.TimeoutExpired:
                raise EngineError(f"{exe} did not finish within {timeout} s") from None
            log += p.stdout + p.stderr
            if p.returncode != 0:
                raise EngineError(f"{exe} exited with status {p.returncode}\n{log[-3000:]}")
    elif runner == "docker":
        image = DOCKER_IMAGE.format(version=version)
        script = " && ".join(" ".join(s) for s in steps)
        try:
            p = subprocess.run(["docker", "run", "--rm", "--platform", "linux/amd64",
                                "-v", f"{work_dir}:/work", "-w", "/work", image,
                                "bash", "-c", script],
                               capture_output=True, text=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            raise EngineError(f"docker run {image} did not finish within {timeout} s") from None
        log = p.stdout + p.stderr
        if p.returncode != 0:
            raise EngineError(f"docker run {image} exited with status {p.returncode}\n{log[-3000:]}")
    else:
        raise ValueError(f"unknown runner {runner!r}")

    err = _engine_error(log)
    if err or _FINISHED not in log:
        raise EngineError(f"HEC-RAS {version} did not finish: {err or 'no completion message'}\n"
                          f"{log[-3000:]}")
    os.replace(work_dir / tmp_hdf, out_hdf)
    m = _VOLUME_ERR.search(log)
    return RunResult(out_hdf, log, float(m.group(1)) if m else None)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("work_dir", type=Path)
    ap.add_argument("project")
    ap.add_argument("--plan", default="01")
    ap.add_argument("--version", default=os.environ.get("HECRAS_VERSION", DEFAULT_VERSION))
    ap.add_argument("--runner", choices=["native", "docker"],
                    default=os.environ.get("HECRAS_RUNNER") or None)
    args = ap.parse_args(argv)
    try:
        res = run_engine(args.work_dir, args.project, args.plan, args.version, args.runner)
    except (EngineError, FileNotFoundError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    print(f"Results: {res.plan_hdf}")
    if res.volume_error_pct is not None:
        print(f"Volume accounting error: {res.volume_error_pct:.5f} %")
    return 0


if __name__ == "__main__":
    sys.exit(main())
