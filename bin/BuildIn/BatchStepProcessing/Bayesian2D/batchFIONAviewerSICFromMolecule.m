function Report = batchFIONAviewerSICFromMolecule(Molecule, varargin)
%BATCHFIONAVIEWERSICFROMMOLECULE Run FIONAviewer batch fitting from Molecule.
%   REPORT = BATCHFIONAVIEWERSICFROMMOLECULE(MOLECULE) builds temporary
%   FIONA-format traces directly from the loaded FIESTA Molecule structures,
%   invokes batchFIONAviewerSIC on them, and writes final results to a folder
%   selected by the user. No pre-existing *_fiona.mat folder is required.

    p = inputParser;
    p.addRequired('Molecule', @isstruct);
    p.addParameter('OutputFolder', '', @(x) ischar(x) || isstring(x));
    p.addParameter('OutputFolderName', 'FIONAviewer_SIC5_FromMolecule', @(x) ischar(x) || isstring(x));
    p.addParameter('OnlySelected', false, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
    p.addParameter('PathDegree', 1, @(x) isnumeric(x) && isscalar(x) && ismember(x,[1 2]));
    p.addParameter('DeleteTemporaryFiles', true, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
    p.KeepUnmatched = true;
    p.parse(Molecule, varargin{:}); opt = p.Results;
    if isempty(opt.OutputFolder)
        outputParent = uigetdir(pwd, 'Select parent folder for processed FIONAviewer traces');
        if isequal(outputParent,0), Report = struct([]); return, end
        outputFolder = fullfile(outputParent, char(opt.OutputFolderName));
    else
        outputFolder = char(opt.OutputFolder);
    end
    if ~isfolder(outputFolder), mkdir(outputFolder); end

    thisFolder = fileparts(mfilename('fullpath'));
    addpath(fileparts(thisFolder), '-begin'); % projection/path helper functions
    indices = 1:numel(Molecule);
    if opt.OnlySelected
        if ~isfield(Molecule,'Selected'), error('batchFIONAviewerSICFromMolecule:Selected', 'Molecule has no Selected field.'); end
        indices = find([Molecule.Selected] == 1);
    end
    temporaryFolder = tempname;
    mkdir(temporaryFolder);
    cleanup = onCleanup(@() cleanupTemporaryFolder(temporaryFolder, opt.DeleteTemporaryFiles));
    for q = 1:numel(indices)
        k = indices(q);
        data = moleculeToFIONAData(Molecule(k), opt.PathDegree);
        save(fullfile(temporaryFolder, [moleculeFileStem(Molecule(k),k) '_fiona.mat']), 'data', '-v7');
    end
    if isempty(indices), Report = struct([]); return, end

    % Forward all batchFIONAviewerSIC options (ShowViewer, StagePauseSeconds,
    % DownsampleFactor, etc.) while forcing its final output folder.
    unmatched = namedArgsToCell(p.Unmatched);
    Report = batchFIONAviewerSIC(temporaryFolder, unmatched{:}, 'OutputFolder', outputFolder);
    clear cleanup
end

function data = moleculeToFIONAData(Molecule, pathDegree)
    if ~isfield(Molecule,'Results') || size(Molecule.Results,2) < 4
        error('batchFIONAviewerSICFromMolecule:Results', 'Molecule lacks Results(:,2:4).');
    end
    results = double(Molecule.Results);
    time = results(:,2); xy = results(:,3:4);
    % Match exportStepFit2DToFIONA: use the 2-D pipeline's stored path when
    % available, otherwise derive a fresh linear/quadratic PCA path from the
    % raw x-y trace. PathData is deliberately not used because its columns
    % have different meanings across FIESTA workflows.
    if isfield(Molecule,'PathFit2D') && ~isempty(Molecule.PathFit2D)
        model = Molecule.PathFit2D;
    else
        model = fitStepFit2DPath(xy,pathDegree);
    end
    [longAxis,shortAxis] = projectStepFit2DPath(xy,model);
    longAxis = longAxis - longAxis(1);
    data.xy = [longAxis shortAxis];
    data.yx = [shortAxis longAxis];
    data.time = time;
    data.neighbors = [];
    if size(results,2) >= 7, data.efo = results(:,7); end
    if isfield(Molecule,'StepFit2D') && isstruct(Molecule.StepFit2D) && ...
            isfield(Molecule.StepFit2D,'EFOTruncation') && isfield(Molecule.StepFit2D.EFOTruncation,'Threshold')
        data.efoThreshold = Molecule.StepFit2D.EFOTruncation.Threshold;
    end
end

function stem = moleculeFileStem(Molecule,index)
    stem = sprintf('%d',index);
    if isfield(Molecule,'Name') && ~isempty(Molecule.Name)
        stem = Molecule.Name;
        if iscell(stem), stem = stem{1}; end
        if isnumeric(stem), stem = num2str(stem); elseif isstring(stem), stem = char(stem); end
        stem = regexprep(strtrim(stem), '^Molecule\s*', '', 'ignorecase');
    end
    stem = regexprep(stem, '[^A-Za-z0-9_-]+', '_');
    if isempty(stem), stem = sprintf('%d',index); end
end

function args = namedArgsToCell(values)
    fields = fieldnames(values); args = cell(1,2*numel(fields));
    for j = 1:numel(fields), args{2*j-1} = fields{j}; args{2*j} = values.(fields{j}); end
end

function cleanupTemporaryFolder(folder, shouldDelete)
    if shouldDelete && isfolder(folder), rmdir(folder,'s'); end
end
