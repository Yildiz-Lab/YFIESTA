% RUNMOLECULEBAYES2DPIPELINE
% Run after loading a MAT file containing Molecule (and optionally Filament).
% This package uses global Bayesian MAP segmentation, not greedy SIC splits.

packageFolder = fileparts(mfilename('fullpath'));
addpath(packageFolder);
addpath(fileparts(packageFolder), '-begin'); % shared helpers; prefer this package's current versions
if ~exist('Molecule','var')
    error('Load a MAT file containing the variable ''Molecule'' before running this script.');
end

Options = bayesStepFit2DDefaultOptions();
Options.Bayes.DecimationFactor = 10;
Options.Bayes.ExpectedDwellPoints = 30; % larger = more conservative
Options.Path.Degree = 1;

[MoleculeBayesian2D,MoleculeBayesianPassed,Bayesian2DReport,Bayesian2DOptions] = ...
    processMoleculeBayes2D(Molecule,Options);
save('Molecule_Bayesian2D.mat','MoleculeBayesian2D','Bayesian2DReport','Bayesian2DOptions','-v7.3');
save('Molecule_Bayesian2D_passed_truncated.mat','MoleculeBayesianPassed', ...
    'Bayesian2DReport','Bayesian2DOptions','-v7.3');

% FIESTA-compatible subset, using native SaveTracks variable names/schema.
Molecule = MoleculeBayesian2D;
passedMask = arrayfun(@(m) isstruct(m.StepFit2DStatus) && ...
    isfield(m.StepFit2DStatus,'Passed') && m.StepFit2DStatus.Passed,Molecule);
Molecule(~passedMask) = [];
if ~exist('Filament','var'), Filament = []; end
saveStepFit2DForFiesta('Molecule_Bayesian2D_FIESTA.mat',Molecule,Filament);

FIONAExportReport = exportStepFit2DToFIONA(MoleculeBayesianPassed,'FIONA_Bayesian2D');
save('Molecule_Bayesian2D_FIONA_export_report.mat','FIONAExportReport','-v7');
fprintf('Saved %d Bayesian-passed molecules out of %d total.\n',numel(MoleculeBayesianPassed),numel(MoleculeBayesian2D));
