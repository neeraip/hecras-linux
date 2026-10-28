# Example projects

Each folder is a HEC-RAS 6.6 project saved by the Windows GUI. All are 2D-only
(one 2D flow area, no structures) with US Customary units.

| Folder | Project | Interior cells | Boundary conditions | Infiltration | Runs as shipped |
|---|---|---|---|---|---|
| `AFTER_RUN` | `test_hdf` | 375 | Flow hydrograph on `US`, normal depth on `DS` | yes | Workflow A |
| `BEFORE_RUN` | `test_hdf` | 375 | Same project before "Compute Geometry": mesh only, no tables, no results, no `.b01`/`.x01` (borrow them and the BC data from `AFTER_RUN`) | | Workflow B |
| `BEC` | `BEC` | 2,315 | Flow hydrograph on `DS`, normal depth on `Inflow` (as named in the project) | SCS curve number | Workflow A |
| `BEC_WO_Infiltration` | `model1` | 2,622 | Flow hydrograph on `inflow`, normal depth on `DS` | no | Workflow A |
| `VA` | `VA` | 1,221 | Normal depth on `BC Line 1` (around the area) + precipitation | no | Workflow A |
| `Muncie` | `Muncie` | 3,146 | Flow hydrograph on `US`, normal depth on `DS` + precipitation | yes | Workflow A |

## Using them for Workflow B

The GUI-computed projects double as test inputs for the other workflows, with
the GUI results as the reference:

```bash
python3 tools/make_workflow_input.py Examples/BEC /tmp/BEC_B --workflow B   # mesh only
python3 ras_preprocess.py /tmp/BEC_B --output-dir /tmp/BEC_B_run \
        --event-conditions Examples/BEC/BEC.p01.hdf
cp Examples/BEC/BEC.b01 Examples/BEC/BEC.x01 /tmp/BEC_B_run/
python3 hecras_run.py /tmp/BEC_B_run BEC
python3 tools/compare_results.py /tmp/BEC_B_run/BEC.p01.hdf Examples/BEC/BEC.p01.hdf
```

## DSS files are not included

The flow hydrographs of these projects live in DSS files, which `.gitignore`
excludes. Workflow A does not need them: the GUI's `.p01.hdf` already holds
the boundary data. Workflow B builds the plan from scratch, so without
the DSS file the run has no inflow. Pass `--event-conditions <GUI .p01.hdf>`
to reuse the GUI's boundary data, as the tests do.
