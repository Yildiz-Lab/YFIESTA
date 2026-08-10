function ExportReport = exportStepFit2DToFIONA(Molecules, outputFolder)
%EXPORTSTEPFIT2DTOFIONA Export passed 2-D fits in FIONAviewer MAT format.
%   exportStepFit2DToFIONA(MOLECULES, FOLDER) writes one file per molecule:
%   <Molecule.Name>_fiona.mat.  The data structure matches FIONAviewer's
%   expected schema and includes both long- and short-axis step fits.
%
%   data.trace columns are:
%       [longRaw shortRaw longFit shortFit longChangePoint shortChangePoint]
%   data.trace_yx keeps the short-axis fit in FIONAviewer's alternate view.

    if nargin < 2 || isempty(outputFolder)
        outputFolder = fullfile(pwd, 'FIONA_2D');
    end
    if ~isfolder(outputFolder)
        mkdir(outputFolder);
    end
    ExportReport = repmat(struct('Name','','File','','Succeeded',false,'Reason',''), numel(Molecules), 1);

    for k = 1:numel(Molecules)
        ExportReport(k).Name = fionaName(Molecules(k), k);
        try
            if ~isfield(Molecules(k), 'StepFit2D') || isempty(Molecules(k).StepFit2D)
                error('Missing StepFit2D result.');
            end
            if ~isfield(Molecules(k), 'PathFit2D') || isempty(Molecules(k).PathFit2D)
                error('Missing PathFit2D model.');
            end
            results = Molecules(k).Results;
            if size(results,2) < 4
                error('Results must contain time, x, and y columns.');
            end
            rawXY = double(results(:,3:4));
            time = double(results(:,2));
            stepFitXY = double(Molecules(k).StepFit2D.StepFit);
            if ~isequal(size(rawXY), size(stepFitXY))
                error('StepFit2D.StepFit does not match the truncated Results length.');
            end

            pathModel = Molecules(k).PathFit2D;
            [longRaw, shortRaw] = projectStepFit2DPath(rawXY, pathModel);
            [longFit, shortFit] = projectStepFit2DPath(stepFitXY, pathModel);

            % FIONA's long-axis convention begins at the first raw point;
            % use that same reference for the long-axis fitted plateaus.
            longReference = longRaw(1);
            longRaw = longRaw - longReference;
            longFit = longFit - longReference;

            longCP = [0; double(diff(longFit) ~= 0)];
            shortCP = [0; double(diff(shortFit) ~= 0)];
            data.xy = [longRaw shortRaw];
            data.yx = [shortRaw longRaw];
            data.neighbors = [];
            data.time = time;
            data.trace = [longRaw shortRaw longFit shortFit longCP shortCP];
            data.trace_yx = nan(numel(time), 6);
            data.trace_yx(:,[1 3 5]) = [shortRaw shortFit shortCP];
            % Preserve EFO and the earlier robust truncation threshold when
            % available, so downstream manual-review tools can flag high-EFO
            % regions without reopening the original Molecule MAT file.
            if size(results,2) >= 7
                data.efo = double(results(:,7));
            end
            if isfield(Molecules(k).StepFit2D, 'EFOTruncation') && ...
                    isfield(Molecules(k).StepFit2D.EFOTruncation, 'Threshold')
                data.efoThreshold = double(Molecules(k).StepFit2D.EFOTruncation.Threshold);
            end
            if isfield(Molecules(k).StepFit2D, 'EFOFlags') && ...
                    isfield(Molecules(k).StepFit2D.EFOFlags, 'HighMask')
                data.efoHigh = logical(Molecules(k).StepFit2D.EFOFlags.HighMask(:));
                data.efoHighRegions = Molecules(k).StepFit2D.EFOFlags.Regions;
            end

            filename = [fionaName(Molecules(k), k) '_fiona.mat'];
            fullFilename = fullfile(outputFolder, filename);
            save(fullFilename, 'data', '-v7');
            ExportReport(k).File = fullFilename;
            ExportReport(k).Succeeded = true;
            ExportReport(k).Reason = 'Exported';
        catch ME
            ExportReport(k).Reason = ME.message;
            warning('exportStepFit2DToFIONA:ExportFailed', ...
                'Could not export %s: %s', ExportReport(k).Name, ME.message);
        end
    end
end

function name = fionaName(Molecule, index)
    name = sprintf('%d', index);
    if isfield(Molecule, 'Name') && ~isempty(Molecule.Name)
        name = char(Molecule.Name);
        name = regexprep(strtrim(name), '^Molecule\s*', '', 'ignorecase');
    end
    % Preserve an integer molecule number when available, but keep a safe
    % filename for other FIESTA Name values.
    name = regexprep(name, '[^A-Za-z0-9_-]+', '_');
    if isempty(name), name = sprintf('%d', index); end
end
