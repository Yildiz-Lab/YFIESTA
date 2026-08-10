# FIONAviewer batch SIC wrapper

> This folder deliberately contains only the folder-based FIONAviewer batch
> workflow. For raw `Molecule` data, use
> `runMolecule2DSICThenFIONAviewerSIC` in the parent `outputs/` folder. It
> first runs the established 2-D SIC preparation/fitting pipeline, exports
> its passed traces, and then invokes this folder-based batch.

Run `runFIONAviewerBatchSIC` in MATLAB to select a folder containing
`*_fiona.mat` files. Alternatively, use
`runFIONAviewerBatchSIC('/path/to/fiona/files')`. The wrapper creates
`FIONAviewer_SIC5/` in that folder.
For every input trace it uses the existing FIONAviewer functions to decimate
by five, fit and align the long axis, then fit and align the short axis.

The output subfolder receives the updated MAT file, an exported `.fig`, and a
`FIONAviewer_SIC5_batch_report.mat` status report. After final FIONA scoring,
each trace and its figures are organized into `S/`, `A/`, `B/`, or
`C_or_lower/`, matching the 2-D SIC organization. `all/`
contains one MAT copy per final trace and is used automatically by the Step
Statistics GUI when its parent FIONA output folder is selected. Set
`ShowViewer` to true when calling `batchFIONAviewerSIC` directly to watch each
trace in the GUI.

## Steppiness review reports

The batch also adds `data.Steppiness` to every output MAT file and writes:

- `steppiness_trace_report.csv`: one A–F score per trace.
- `steppiness_region_report.csv`: one score per fitted dwell region.
- `<trace>_review_flags.fig`: long-/short-axis review plot; black is
  unflagged, red is a weak/short-step region, and blue is high EFO.

Trace grades run from `S` (best) through `A`, `B`, `C`, `D`, and `F`. The
heuristic combines step evidence, trace length, fitted-path directionality,
and how closely the fitted long-axis progression follows constant velocity.
It separately flags short/weak step regions and high-EFO regions. These are
review priorities rather than automatic biological labels.

High-EFO flags require `data.efo` and `data.efoThreshold`, which are added to
new FIONA files exported by the 2-D pipeline. Older FIONA files remain fully
gradable for steppiness, but their report will show `EFOAvailable = false`.

## Optional parameters

`batchFIONAviewerSIC` accepts these name/value options:

| Parameter | Default | Meaning |
|---|---:|---|
| `DownsampleFactor` | `5` | FIONAviewer **Filter Data → Decimate → Dec/Window** value used before fitting each axis. |
| `OutputFolderName` | `'FIONAviewer_SIC5'` | Name of the new results subfolder created inside the selected input folder. |
| `ShowViewer` | `false` | Set `true` to show FIONAviewer while each long-/short-axis filtering, fitting, and alignment stage is performed. |
| `StagePauseSeconds` | `1` when `ShowViewer=true`; otherwise `0` | Time to pause at every visible processing stage. |

Examples:

```matlab
% Select the folder interactively and watch each processing stage.
batchFIONAviewerSIC('', 'ShowViewer', true, 'StagePauseSeconds', 2)

% Process a known folder without opening the GUI, using 10x decimation.
batchFIONAviewerSIC('/path/to/fiona/files', 'DownsampleFactor', 10)

% Use the convenience wrapper, which opens uigetdir when no folder is given.
runFIONAviewerBatchSIC([], 'ShowViewer', true)
```
