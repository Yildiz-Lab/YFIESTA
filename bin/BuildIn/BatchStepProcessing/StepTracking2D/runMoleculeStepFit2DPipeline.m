% RUNMOLECULESTEPFIT2DPIPELINE
% Run this script after loading a MAT file containing the structure array
% named Molecule.  Edit the options below during validation.

if ~exist('Molecule', 'var')
    error('Load a MAT file containing the variable ''Molecule'' before running this script.');
end

Options = stepFit2DDefaultOptions();
Options.Fit.DecimationFactor = 10;
Options.Path.Degree = 1;             % Set to 2 for a curved filament/path.
Options.Plot.InitialFits = false;
Options.Plot.EFOAfterTruncation = true;
Options.Plot.AlignedFits = true;

[MoleculeStepFit2D, MoleculePassedTruncated, StepFit2DReport, StepFit2DOptions] = ...
    processMoleculeStepFit2D(Molecule, Options);

save('Molecule_StepFit2D.mat', 'MoleculeStepFit2D', 'StepFit2DReport', 'StepFit2DOptions', '-v7.3');
save('Molecule_StepFit2D_passed_truncated.mat', 'MoleculePassedTruncated', ...
    'StepFit2DReport', 'StepFit2DOptions', '-v7.3');

% GUI-compatible output: mirror FIESTA SaveTracks by starting with the
% complete structure array and deleting all non-passing molecules in place.
% This preserves the original FIESTA struct-array orientation and schema.
Molecule = MoleculeStepFit2D;
passedMask = arrayfun(@(m) isstruct(m.StepFit2DStatus) && ...
    isfield(m.StepFit2DStatus, 'Passed') && m.StepFit2DStatus.Passed, Molecule);
Molecule(~passedMask) = [];
if ~exist('Filament', 'var')
    Filament = [];
end
saveStepFit2DForFiesta('Molecule_StepFit2D_FIESTA.mat', Molecule, Filament);

% One FIONAviewer-compatible MAT file per passed molecule, including both
% longitudinal and transverse projections and their 2-D-derived step fits.
FIONAExportReport = exportStepFit2DToFIONA(MoleculePassedTruncated, 'FIONA_2D');
save('Molecule_StepFit2D_FIONA_export_report.mat', 'FIONAExportReport', '-v7');

fprintf('Saved %d passed molecules out of %d total.\n', ...
    numel(MoleculePassedTruncated), numel(Molecule));
