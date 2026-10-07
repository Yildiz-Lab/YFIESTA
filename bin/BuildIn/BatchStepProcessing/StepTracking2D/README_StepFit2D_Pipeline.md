# 2-D Molecular-Motor Step-Fitting Pipeline

Run `runMoleculeStepFit2DPipeline.m` after loading a MAT file that contains
the FIESTA structure array `Molecule`.

The wrapper performs this sequence for every molecule:

1. Quality filtering (`assessStepFit2DMolecule`) using finite coordinates,
   track length, net displacement, directionality, and PCA elongation.
2. Initial EFO truncation (`truncateHighEFOStepFit2D`): it retains the first
   run of five points at or below the global median EFO (a robust baseline
   estimate).
3. Polynomial path modelling (`fitStepFit2DPath`): set `Options.Path.Degree`
   to 1 for a straight path or 2 for a curved path.
4. Path-constrained SIC segmentation (`SICstepFinder2D`) on binned/decimated
   data. The default uses the fitted path's longitudinal coordinate to choose
   changepoints, with global dynamic programming to find the lowest-SIC
   segmentation across all possible numbers of dwells. Full x-y dwell means
   are then refit and retained in the output.
5. Initial plotting, longitudinal local-displacement alignment
   (`alignStepFit2D`), and replotting. A final raw longitudinal SIC comparison
   removes an aligned boundary only when merging it improves the model. Thus,
   short transverse blips are not removed by a hard lifetime threshold.
6. Saving of all input molecules with status fields, plus the passed and
   EFO-truncated subset.

## Saved files

- `Molecule_StepFit2D.mat`: `MoleculeStepFit2D` includes every original
  molecule. Passed molecules have truncated `Results`, `StepFit2D`, and
  `PathFit2D`; rejected molecules retain raw data and a `StepFit2DStatus`.
- `Molecule_StepFit2D_passed_truncated.mat`: `MoleculePassedTruncated` holds
  only passed, EFO-truncated molecules.
- `Molecule_StepFit2D_FIESTA.mat`: GUI-compatible output with variables
  named exactly `Molecule` and `Filament`. It contains only quality/EFO-
  passed, truncated molecules. Load this file through FIESTA's **Load
  Tracks** command for visualization. It omits pipeline-only fit fields to
  remain compatible with FIESTA's structure concatenation.
- `FIONA_2D/`: one FIONAviewer-compatible `<number>_fiona.mat` file per
  passed molecule. Each includes long-/short-axis raw projections, both
  stepwise fits, and changepoints in the `data.trace` convention expected
  by FIONAviewer.

When using `batchMolecule2DSICThenFIONAviewerSIC`, the corresponding
`FIONA_Input_2DSIC/` folder is organized after 2-D scoring into `S/`, `A/`,
`B/`, and `C_or_lower/`. These grades use the same
`scoreFIONASteppiness` evaluation as the FIONAviewer stage, applied to the
2-D projected fit. Each grade folder contains the corresponding
`*_fiona.mat`, `*_2DSIC_trajectory.fig`, and `*_2DSIC_projection.fig` files.
`For_FIONAviewer_SIC/` is a staging folder containing
only traces at or above the configured pre-screen score (default 40, i.e.
S through D), so the normal `FIONAviewer_SIC/` batch continues to run while
F-grade traces remain available for manual review in `C_or_lower/`.

## SIC segmentation options

The default `Options.Fit.SegmentationMode = 'pathConstrained'` is intended
for processive motors. It suppresses changepoints driven solely by off-axis
noise and avoids the greedy algorithm's tendency to omit two nearby on-axis
steps. Set it to `'joint2D'` to reproduce the earlier greedy 2-D SIC
behavior:

```matlab
Options = stepFit2DDefaultOptions;
Options.Fit.SegmentationMode = 'joint2D';
```

`DecimationFactor = 10` averages ten consecutive localizations only for the
initial SIC search; fitted values are expanded and plotted over the original
data. The final alignment returns to the raw data and shifts each boundary
within plus/minus `DecimationFactor` samples to the maximal longitudinal
mean displacement, analogous to FIONAviewer's `AlignToMaxDx`. The default
`Alignment.SearchWindow = []` enables that automatic coupling; set a numeric
value to override it. `MinPointsInStep = 2` is retained as the numerical
minimum needed to estimate a variance from a decimated dwell. It is not a
biological lifetime cutoff, and the post-alignment pruning uses SIC evidence
rather than a new minimum-duration rule.

### The model does not encode a 16 nm step size

No default option specifies 16 nm, a fixed step size, or an integer step
lattice. Dwell positions are estimated from the data, and SIC selects their
number and boundaries. Consequently, approximately 16 nm spacing can emerge
when supported by the longitudinal data, but it is neither forced nor used to
reject a differently sized step. This is deliberate: motors, calibration,
and path projection can all shift apparent step sizes.

## Tune before interpreting results

The automatic filter is deliberately exposed in `stepFit2DDefaultOptions`.
Most importantly, validate and adjust `MinPoints`, `MinNetDisplacement`,
`MinDirectionality`, `MinElongation`, the EFO rule, `DecimationFactor`, and
`SegmentationMode` against manually curated trajectories before using it as
an exclusion rule. The default `MinDirectionality` is `0.002`, selected to
retain slow, elongated traces whose long dwells otherwise inflate the
point-to-point path length.

For quadratic paths, nearest-path projection is more computationally costly
than a linear path because each localization is projected to the curve.

## Analyze a collection of FIESTA files

Use `batchFiestaCollection2DSICThenFIONAviewerSIC` when selected FIESTA MAT
files from several acquisition folders should be processed and reviewed as
one experiment. Pass a cell array of full file paths and one collection
output folder:

```matlab
[CollectionReport, CombinedManifest] = ...
    batchFiestaCollection2DSICThenFIONAviewerSIC(sourceFiles, outputFolder);
```

The wrapper copies the selected inputs into `Source_FIESTA/`, processes each
file independently below `Per_Dataset/`, and gathers the final outputs below
`Combined_FIONAviewer_SIC/S`, `A`, `B`, `C_or_lower`, and `all`. Combined
filenames are prefixed with their source-dataset name so identical molecule
numbers from different acquisitions do not collide. Progress is written
after every dataset, and the default `ResumeExisting=true` reuses datasets
that have a `collection_dataset_complete.mat` marker.

## Export S/A step sizes and dwell times

`exportFIONASAStatsCSV` collects the final aligned FIONAviewer results from
the `S` and `A` grade folders. It writes two columns (`step size (nm)` and
`dwell times (ms)`) and leaves one blank row between molecules. Pass the
original overarching FIESTA MAT file, any enclosing analysis folder, or the
`FIONAviewer_SIC` folder itself:

```matlab
Report = exportFIONASAStatsCSV('/full/path/to/dataset_fiesta.mat');
```

The default output is
`FIONAviewer_SIC/steps_summary_S_A.csv`. Call the function without an input
to select the overarching MAT file with a file-selection dialog, or provide
a second argument to choose another CSV filename.
