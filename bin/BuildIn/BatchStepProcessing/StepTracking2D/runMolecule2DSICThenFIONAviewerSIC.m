% RUNMOLECULE2DSICTHENFIONAVIEWERSIC
% Run after loading a MAT file containing the FIESTA variable Molecule.
% Performs 2-D SIC preparation/fitting, exports the passed traces, then
% runs the folder-based FIONAviewer SIC batch and its review evaluation.

if ~exist('Molecule', 'var')
    error('Load a MAT file containing the variable ''Molecule'' before running this script.');
end

if exist('Filament', 'var')
    batchMolecule2DSICThenFIONAviewerSIC(Molecule, 'Filament', Filament);
else
    batchMolecule2DSICThenFIONAviewerSIC(Molecule);
end
