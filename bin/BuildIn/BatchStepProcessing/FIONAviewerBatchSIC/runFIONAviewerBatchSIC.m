function Report = runFIONAviewerBatchSIC(varargin)
%RUNFIONAVIEWERBATCHSIC Run FIONAviewer batch fitting on a selected folder.
%   runFIONAviewerBatchSIC() opens a folder-selection dialog.
%   runFIONAviewerBatchSIC(FOLDER) processes FOLDER directly.

    % Several legacy copies of batchFIONAviewerSIC can coexist on a
    % YFIESTA MATLAB path.  Always resolve the batch implementation next
    % to this runner, rather than accidentally calling an older variant.
    thisFolder = fileparts(mfilename('fullpath'));
    addpath(thisFolder, '-begin');
    if nargin < 1 || isempty(varargin{1})
        inputFolder = uigetdir(pwd, 'Select folder containing *_fiona.mat files');
        if isequal(inputFolder, 0)
            Report = struct([]);
            return
        end
    else
        inputFolder = varargin{1};
    end
    Report = batchFIONAviewerSIC(inputFolder, varargin{2:end});
end
