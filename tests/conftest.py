"""pytest options and engine discovery for the pipeline tests.

    pytest                               # every installed engine version
    pytest --hecras-version 7.0.1        # one version
    pytest --runner docker               # use the hecras-linux:<version> images
    pytest -m "not engine"               # preprocessing only, no HEC-RAS needed
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

import hecras_run  # noqa: E402


def pytest_addoption(parser):
    parser.addoption("--hecras-version", action="append", default=None,
                     help="HEC-RAS engine version(s) to test (default: all under bin/)")
    parser.addoption("--runner", choices=["native", "docker"], default=None,
                     help="how to run the engines (default: native on Linux x86_64)")


def pytest_configure(config):
    config.addinivalue_line("markers", "engine: needs the HEC-RAS Linux engines")
    config.addinivalue_line("markers", "experimental: Workflow B/C, not yet validated")


def _engine_params(config):
    runner = config.getoption("--runner") or hecras_run.default_runner()
    wanted = config.getoption("--hecras-version") or hecras_run.available_versions()
    params = []
    for v in wanted:
        if runner == "native":
            ok = hecras_run.native_supported() and v in hecras_run.available_versions()
            why = f"HEC-RAS {v} engines not runnable natively here"
        else:
            ok = hecras_run.docker_image_present(v)
            why = f"Docker image {hecras_run.DOCKER_IMAGE.format(version=v)} not built"
        marks = [] if ok else [pytest.mark.skip(reason=why)]
        params.append(pytest.param((v, runner), id=f"ras{v}", marks=marks))
    if not params:
        params.append(pytest.param(None, id="no-engine",
                                   marks=pytest.mark.skip(reason="no engines under bin/")))
    return params


def pytest_generate_tests(metafunc):
    if "engine" in metafunc.fixturenames:
        metafunc.parametrize("engine", _engine_params(metafunc.config))
