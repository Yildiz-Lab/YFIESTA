
function fConvertMINFLUX(minflux_file, savefull)

% JS Repackage into Molecule-like thing for FIESTA analysis base


if nargin < 2
    % Now package to save in a separate MAT file
    [savepath,savename,~] = fileparts(minflux_file);
    savefull = fullfile(savepath,strcat(savename,'_fiesta','.mat'));
end


[path, fname, ~] = fileparts(minflux_file);
MFData = load(minflux_file);

% MF data looks is just another structure.

% itr - contains most of the data, some of the important ones we will
% incorporate
%       loc - contains non-drift corrected localization tuples
%       lnc - contains drift corrected localization tuples
%       eco - number of photon counts collected for the (?duration of the
%       iteration scheme OR just the final?)
%       cfr - center frequency ratio, if checked tells how well centered
%       your fluorophore is

% gri
% tim - time of acquisition. Since MINFLUX tracks a single particle at a
    % time, this contains the index of time associated with that data point
% tid - trace ID. A consecutive trace has a specific ID (kind of like a
    % molecule in FIESTA)
% vld - channel...maybe?
% act - ?
% dos - ?
% sky - phase of localization achieved?

% load an example file (because hell it doesn't work if I initialize it
% myself)

example_struct = load('example_molformat.mat');
MFMolecule = example_struct.Molecule;

% First, let's initialize a Molecule struct
uid = unique(MFData.tid);

% MFMolecule = struct();

Colorwheel = [0,0,1; 0,1,0; 1,0,0];

for j = 1:length(uid)

    idx = find(MFData.tid == uid(j));

    MFMolecule(j).Name = ['Molecule ', ' ', num2str(uid(j))];
    MFMolecule(j).File = fname;
    MFMolecule(j).Comments = ['MINFLUX Acq [efo=fwhm, eco=Amp]']; % Comments can also be applied as "MINFLUX Acquisition"
    

    % Molecule.Results is the most important change
    % [Frame, Time, XPos (pix), YPos (pix), ZPos (pix), Distance, FWHM,
    % Amplitude (counts), Position Error (nm), Tags]
    % Results packaging
    MFMolecule(j).Results = nan(length(idx),10);
    MFMolecule(j).Results(:,1) = idx';
    MFMolecule(j).Results(:,2) = MFData.tim(idx)'; %s
    MFMolecule(j).Results(:,3) = MFData.itr.loc(idx,4,1)*10^9; %converted to nm
    MFMolecule(j).Results(:,4) = MFData.itr.loc(idx,4,2)*10^9;
    % MFMolecule(j).Results(:,5) = MFData.itr.loc(:,4,3)*10^9;
    MFMolecule(j).Results(:,5) = nan(length(idx),1);
    MFMolecule(j).Results(:,6) = pdist2(MFMolecule(j).Results(:,3:4), MFMolecule(j).Results(1,3:4));
    MFMolecule(j).Results(:,7) = sum(MFData.itr.efo(idx,:),2,"omitnan"); %this is now the efo
    MFMolecule(j).Results(:,8) = sum(MFData.itr.eco(idx,:),2); %this is now the eco
    MFMolecule(j).Results(:,9) = nan(length(idx),1);
    MFMolecule(j).Results(:,10) = zeros(length(idx),1);


    % Mostly fixed
    MFMolecule(j).Type = 'symmetric';
    MFMolecule(j).Selected = 0;
    MFMolecule(j).Visible = 1;
    MFMolecule(j).Drift = 0;

    % MINFLUX we will set the pixel value to 1 nm, this way there is no
    % conversion necessary for the file. It should also then be another sign
    % you are dealing with a MINFLUX data set
    MFMolecule(j).PixelSize = 1;


    MFMolecule(j).Channel = MFData.vld(idx(1));
    MFMolecule(j).Color = Colorwheel(MFMolecule(j).Channel,:);
    MFMolecule(j).PathData = [];
    MFMolecule(j).PlotHandles = [];
    MFMolecule(j).TrackingResults = {};

end

% clear the rest

for j = length(MFMolecule):-1:length(uid)+1
    MFMolecule(j) = [];
end

Molecule = MFMolecule;

% Make struct field cause it needs it
Filament = [];
Filament=fDefStructure(Filament,'Filament');
Molecule = fDefStructure(Molecule,'Molecule');


% save savefull Config MFMolecule
save(savefull,'Molecule','Filament','-v6');


% For MINFLUX data explicitly, we have a bunch of other parameters that
% might help us debug, such as eco, cfr, etc. These can go after the 10th
% spot in results and be called upon to distinguish MINFLUX vs FIONA traces


end
