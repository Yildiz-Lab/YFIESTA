% RUNMOLECULE2DSICTHENFIONAVIEWERSIC
% Run after loading a MAT file containing the FIESTA variable Molecule.
% Performs 2-D SIC preparation/fitting, exports the passed traces, then
% runs the folder-based FIONAviewer SIC batch and its review evaluation.
% Optional workspace FIONAOptions fields:
%   PopulationEFOGradeCap = 'A' (default), 'B', or 'off'.
%   PopulationEFOMultiplier = 1.5 (default).
%   PopulationEFOReference = [] (eligible median) or calibrated EFO in Hz.

if ~exist('Molecule', 'var')
    error('Load a MAT file containing the variable ''Molecule'' before running this script.');
end

if exist('FIONAOptions', 'var')
    if exist('Filament', 'var')
        batchMolecule2DSICThenFIONAviewerSIC(Molecule, 'Filament', Filament, 'FIONAOptions', FIONAOptions);
    else
        batchMolecule2DSICThenFIONAviewerSIC(Molecule, 'FIONAOptions', FIONAOptions);
    end
elseif exist('Filament', 'var')
    batchMolecule2DSICThenFIONAviewerSIC(Molecule, 'Filament', Filament);
else
    batchMolecule2DSICThenFIONAviewerSIC(Molecule);
end
