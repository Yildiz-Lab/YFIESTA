function Report = batchFIONAviewerSIC(inputFolder, varargin)
%BATCHFIONAVIEWERSIC Batch FIONAviewer SIC fits using its existing routines.
%   REPORT = BATCHFIONAVIEWERSIC() prompts for a folder of *_fiona.mat files.
%   It creates FIONAviewer_SIC5 in that folder, then for each trace:
%     1) loads it through FIONAviewer,
%     2) uses FIONAviewer's Decimate filter (factor 5),
%     3) runs its SICstepFinder on the long axis,
%     4) restores the raw trace and calls its AlignToMaxDx,
%     5) repeats for the short axis,
%     6) saves the MAT file and calls FIONAviewer's ExportFigure.
%
%   This is deliberately orchestration only: filtering, fitting, alignment,
%   trace packaging, plotting, and figure export are all existing
%   FIONAviewer functions.

    p = inputParser;
    p.addOptional('inputFolder', '', @(x) ischar(x) || isstring(x));
    p.addParameter('DownsampleFactor', 5, @(x) isnumeric(x) && isscalar(x) && x >= 1 && x == floor(x));
    p.addParameter('OutputFolderName', 'FIONAviewer_SIC5', @(x) ischar(x) || isstring(x));
    p.addParameter('OutputFolder', '', @(x) ischar(x) || isstring(x));
    p.addParameter('ShowViewer', false, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
    p.addParameter('StagePauseSeconds', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 0));
    p.parse(inputFolder, varargin{:}); opt = p.Results;
    if isempty(opt.StagePauseSeconds)
        if opt.ShowViewer, opt.StagePauseSeconds = 1; else, opt.StagePauseSeconds = 0; end
    end

    if isempty(opt.inputFolder)
        inputFolder = uigetdir(pwd, 'Select folder containing *_fiona.mat files');
        if isequal(inputFolder, 0)
            Report = struct([]); return
        end
    else
        inputFolder = char(opt.inputFolder);
    end
    if ~isfolder(inputFolder)
        error('batchFIONAviewerSIC:Folder', 'Folder does not exist: %s', inputFolder);
    end
    if isempty(opt.OutputFolder)
        outputFolder = fullfile(inputFolder, char(opt.OutputFolderName));
    else
        outputFolder = char(opt.OutputFolder);
    end
    if ~isfolder(outputFolder), mkdir(outputFolder); end

    thisFolder = fileparts(mfilename('fullpath'));
    % Prefer an already configured YFIESTA installation. The fallback finds
    % the enclosing Documents folder, which supports this project's standard
    % Documents/GitHub/YFIESTA layout without relying on the current folder.
    fionaFolder = whichFolder('FIONAviewer.m');
    if isempty(fionaFolder)
        documentsFolder = enclosingFolderNamed(thisFolder, 'Documents');
        if ~isempty(documentsFolder)
            candidate = fullfile(documentsFolder, 'GitHub', 'YFIESTA', 'bin', 'BuildIn', 'FIONAviewer');
            if isfolder(candidate), fionaFolder = candidate; end
        end
    end
    if isempty(fionaFolder)
        error('batchFIONAviewerSIC:FIONAviewer', ...
            'Add the FIONAviewer folder to the MATLAB path, then run again.');
    end
    addpath(fionaFolder, fullfile(fionaFolder, 'SubFns'), '-begin');

    files = dir(fullfile(inputFolder, '*_fiona.mat'));
    % macOS can place AppleDouble sidecar files (for example,
    % ._134_fiona.mat) beside valid MAT files on external drives.  Their
    % names match the normal glob but they are not MATLAB files.
    hiddenFile = startsWith({files.name}, '.');
    ignoredCount = sum(hiddenFile);
    files = files(~hiddenFile);
    if ignoredCount > 0
        fprintf('Ignoring %d hidden macOS metadata file(s) in %s.\n', ...
            ignoredCount, inputFolder);
    end
    [~, order] = sort({files.name}); files = files(order);
    Report = repmat(struct('InputFile','','OutputFile','','Succeeded',false,'Reason','', ...
        'LongAudit',[],'ShortAudit',[],'Steppiness',[],'FlagFigure','', ...
        'Grade','','GradeFolder','','StepStatisticsFile',''),numel(files),1);
    RegionScores = cell(numel(files),1);
    for k = 1:numel(files)
        inputFile = fullfile(inputFolder, files(k).name);
        outputFile = fullfile(outputFolder, files(k).name);
        Report(k).InputFile = inputFile; Report(k).OutputFile = outputFile;
        app = [];
        try
            % Use the existing FIONAviewer loader. The supplied version of
            % FIONAviewer leaves FileName/FilePath uncreated, although its
            % LoadMatDataFile helper expects them; initialize those labels
            % before invoking that unmodified loader.
            app = openFIONATrace(inputFile, opt.ShowViewer);

            % Steps 3 and 5: explicitly reproduce the FIONAviewer UI's
            % Filter Data -> Decimate -> value 5 operation before fitting.
            app = applyFIONADecimate(app, opt.DownsampleFactor);
            [app,longAudit] = runAxisFitAndAlign(app, 'Long-axis', opt.StagePauseSeconds);
            app = applyFIONADecimate(app, opt.DownsampleFactor);
            [app,shortAudit] = runAxisFitAndAlign(app, 'Short-axis', opt.StagePauseSeconds);

            % Existing FIONAviewer export uses StepsFilename to choose both
            % the MAT-file path and its companion .fig output path.
            app.LeftPanelFields.StepsFilename.Text = outputFile;
            app = updateMainPlot(app);
            data = app.data;
            data.BatchFIONAviewerSIC.DownsampleFactor = opt.DownsampleFactor;
            data.BatchFIONAviewerSIC.LongAxis = longAudit;
            data.BatchFIONAviewerSIC.ShortAxis = shortAudit;
            data.Steppiness = scoreFIONASteppiness(data, 'DownsampleFactor', opt.DownsampleFactor);
            save(outputFile, 'data', '-v7');
            ExportFigure(app);
            Report(k).Succeeded = true; Report(k).Reason = 'Saved and exported';
            Report(k).LongAudit = longAudit; Report(k).ShortAudit = shortAudit;
            Report(k).Steppiness = data.Steppiness.Trace;
            try
                [~, baseName] = fileparts(outputFile);
                flagFigure = plotFIONASteppinessFlags(data, 'Visible', 'off', ...
                    'Name', [baseName ' — review flags']);
                flagFigureFile = fullfile(outputFolder, [baseName '_review_flags.fig']);
                % Save while visible: a hidden UIFIGURE/figure can be
                % serialized with Visible='off', making the saved FIG seem
                % unopenable in MATLAB. This mirrors the 2-D SIC export.
                flagFigure.Visible = 'on';
                drawnow;
                saveas(flagFigure, flagFigureFile, 'fig');
                flagFigure.Visible = 'off';
                close(flagFigure);
                Report(k).FlagFigure = flagFigureFile;
            catch flagME
                warning('batchFIONAviewerSIC:FlagPlot', ...
                    'Could not create review plot for %s: %s', files(k).name, flagME.message);
            end
            regions = data.Steppiness.Regions;
            for r = 1:numel(regions), regions(r).OutputFile = outputFile; end
            RegionScores{k} = regions;
            deleteFIONAViewer(app);
        catch ME
            Report(k).Reason = ME.message;
            warning('batchFIONAviewerSIC:TraceFailed', 'Failed for %s: %s', files(k).name, ME.message);
            deleteFIONAViewer(app);
        end
    end
    [Report, RegionScores] = organizeFIONAOutputsByGrade(Report, RegionScores, outputFolder);
    save(fullfile(outputFolder, 'FIONAviewer_SIC5_batch_report.mat'), 'Report', 'opt', '-v7');
    writeFIONASteppinessReports(outputFolder, Report, RegionScores);
end

function [Report, RegionScores] = organizeFIONAOutputsByGrade(Report, RegionScores, outputFolder)
    % Keep final FIONA SIC results organized by their post-fit grade. A
    % MAT-only staging folder supports Step Statistics without accidentally
    % including the batch report or duplicated grade-folder contents.
    statisticsFolder = fullfile(outputFolder, 'all');
    if ~isfolder(statisticsFolder), mkdir(statisticsFolder); end
    for k = 1:numel(Report)
        if ~Report(k).Succeeded, continue, end
        grade = finalFIONAGrade(Report(k));
        if isempty(grade)
            warning('batchFIONAviewerSIC:GradeOrganize', ...
                'Could not determine a final grade for %s; leaving it in the root output folder.', ...
                Report(k).OutputFile);
            continue
        end
        gradeFolder = fullfile(outputFolder, finalGradeFolderName(grade));
        if ~isfolder(gradeFolder), mkdir(gradeFolder); end
        oldOutputFile = Report(k).OutputFile;
        [newOutputFile, moved, message] = moveFIONAOutputFile(oldOutputFile, gradeFolder);
        if ~moved
            warning('batchFIONAviewerSIC:GradeOrganize', ...
                'Could not move %s into %s: %s', oldOutputFile, gradeFolder, message);
            continue
        end
        [~, baseName] = fileparts(oldOutputFile);
        moveFIONAOutputFile(fullfile(outputFolder, [baseName '.fig']), gradeFolder);
        if ~isempty(Report(k).FlagFigure)
            [newFlagFigure, flagMoved] = moveFIONAOutputFile(Report(k).FlagFigure, gradeFolder);
            if flagMoved, Report(k).FlagFigure = newFlagFigure; end
        end
        Report(k).OutputFile = newOutputFile;
        Report(k).Grade = grade;
        Report(k).GradeFolder = gradeFolder;
        Report(k).StepStatisticsFile = fullfile(statisticsFolder, [baseName '.mat']);
        [copied, copyMessage] = copyfile(newOutputFile, Report(k).StepStatisticsFile, 'f');
        if ~copied
            warning('batchFIONAviewerSIC:StepStatisticsCopy', ...
                'Could not stage %s for Step Statistics: %s', newOutputFile, copyMessage);
            Report(k).StepStatisticsFile = '';
        end
        if k <= numel(RegionScores) && ~isempty(RegionScores{k})
            for r = 1:numel(RegionScores{k})
                RegionScores{k}(r).OutputFile = newOutputFile;
            end
        end
    end
end

function grade = finalFIONAGrade(reportEntry)
    grade = '';
    if ~isstruct(reportEntry.Steppiness) || ~isfield(reportEntry.Steppiness, 'Grade')
        return
    end
    grade = upper(char(reportEntry.Steppiness.Grade));
    if ~ismember(grade, {'S','A','B','C','D','F'}), grade = ''; end
end

function folder = finalGradeFolderName(grade)
    switch upper(char(grade))
        case 'S', folder = 'S';
        case 'A', folder = 'A';
        case 'B', folder = 'B';
        otherwise, folder = 'C_or_lower';
    end
end

function [destination, moved, message] = moveFIONAOutputFile(source, destinationFolder)
    destination = source;
    moved = false;
    message = '';
    if isempty(source) || ~isfile(source)
        message = 'source file does not exist';
        return
    end
    [~, name, extension] = fileparts(source);
    destination = fullfile(destinationFolder, [name extension]);
    if strcmp(source, destination)
        moved = true;
        return
    end
    [moved, message] = movefile(source, destination, 'f');
end

function [app,audit] = runAxisFitAndAlign(app, axisName, stagePause)
    % Invoke the actual FIONAviewer Fit and Align button callbacks. This is
    % intentionally not a reimplementation of their private onFit/onAlign
    % routines; it is the same sequence as clicking the GUI buttons.
    app.FitPanelFields.DropFitAxis.Value = axisName;
    showBatchStage(app, [axisName ' — decimated trace'], stagePause);
    if strcmp(axisName, 'Long-axis')
        lineObj = findobj(app.AxMain, 'Type', 'Line', 'Tag', 'hLongAxis');
    else
        lineObj = findobj(app.AxMain, 'Type', 'Line', 'Tag', 'hShortAxis');
    end
    if numel(lineObj) ~= 1
        error('batchFIONAviewerSIC:AxisLine', 'FIONAviewer did not create exactly one %s line.', axisName);
    end
    audit.Axis = axisName;
    audit.DecimatedPointCount = numel(lineObj.YData);
    app = pressFIONAButton(app, app.FitPanelFields.BtnFit);
    if strcmp(axisName,'Long-axis')
        audit.StepsBeforeAlignment = sum(diff(app.Data.stepVector) ~= 0);
    else
        audit.StepsBeforeAlignment = sum(diff(app.Data.shortstepVector) ~= 0);
    end
    showBatchStage(app, [axisName ' — SIC fit on decimated trace'], stagePause);

    % Toggle the existing Filter button off, as in the normal UI workflow,
    % so AlignToMaxDx acts on the original trace with FilterWindow still 5.
    app = turnOffFIONAFilter(app);
    showBatchStage(app, [axisName ' — original trace before alignment'], stagePause);
    app = pressFIONAButton(app, app.FitPanelFields.BtnAutoAlign);
    if strcmp(axisName,'Long-axis')
        audit.StepsAfterAlignment = sum(diff(app.Data.stepVector) ~= 0);
    else
        audit.StepsAfterAlignment = sum(diff(app.Data.shortstepVector) ~= 0);
    end
    showBatchStage(app, [axisName ' — aligned original trace'], stagePause);
end

function app = applyFIONADecimate(app, factor)
    % Exact batch equivalent of selecting the Decimate radio button, entering
    % the requested factor, and pressing FIONAviewer's Filter button.
    if app.FilterPanelFields.BtnFilter.UserData.pressed
        app = turnOffFIONAFilter(app);
    end
    app.FilterPanelFields.FilterWindow.Value = factor;
    decimateButton = findobj(app.FilterPanelFields.FilterGroup, 'Text', 'Decimate');
    if numel(decimateButton) ~= 1
        error('batchFIONAviewerSIC:DecimateControl', 'Could not find FIONAviewer''s Decimate control.');
    end
    app.FilterPanelFields.FilterGroup.SelectedObject = decimateButton;
    app = pressFIONAButton(app, app.FilterPanelFields.BtnFilter);
    if ~app.Data.FilteredFlag || numel(app.Data.PSD1Data_Long_Filt) ~= numel(1:factor:numel(app.Data.PSD1Data_Long))
        error('batchFIONAviewerSIC:DecimationFailed', 'FIONAviewer decimation did not produce the expected factor-%d trace.', factor);
    end
end

function app = turnOffFIONAFilter(app)
    if app.FilterPanelFields.BtnFilter.UserData.pressed
        app = pressFIONAButton(app, app.FilterPanelFields.BtnFilter);
    end
end

function app = pressFIONAButton(app, button)
    % Execute the exact callback attached by FIONAviewer.m and retrieve the
    % newly stored application state, just as an interactive click does.
    setappdata(app.fig, 'app', app);
    callback = button.ButtonPushedFcn;
    callback(button, []);
    app = getappdata(app.fig, 'app');
end

function showBatchStage(app, label, seconds)
    if strcmp(app.fig.Visible, 'on')
        app.FileName.Text = label;
        title(app.AxMain, label, 'Interpreter', 'none');
        drawnow;
        if seconds > 0, pause(seconds); end
    end
end

function folder = whichFolder(filename)
    pathToFile = which(filename);
    if isempty(pathToFile), folder = ''; else, folder = fileparts(pathToFile); end
end

function app = openFIONATrace(inputFile, showViewer)
    app = FIONAviewer();
    try
        % Hide immediately, before the loader and callbacks run, so headless
        % batches do not briefly show one FIONAviewer window per trace.
        if ~showViewer, app.fig.Visible = 'off'; end
        app = ensureFIONAFileLabels(app);
        [inputFolder, inputName, inputExtension] = fileparts(inputFile);
        app = LoadMatDataFile(app, inputFolder, [inputName inputExtension]);
        app = updateMainPlot(app);
        setappdata(app.fig, 'app', app);
    catch ME
        % If a loader error occurs before this helper returns, its caller
        % never receives app and therefore cannot close this new UIFIGURE.
        % Dispose of it here before reporting the original failure.
        deleteFIONAViewer(app);
        rethrow(ME)
    end
end

function deleteFIONAViewer(app)
    % UIFIGURE instances are most reliably disposed of with delete rather
    % than close. This guarantees one trace cannot leave a GUI behind.
    if isempty(app) || ~isstruct(app) || ~isfield(app,'fig') || isempty(app.fig)
        return
    end
    if isgraphics(app.fig)
        delete(app.fig);
    end
end

function app = ensureFIONAFileLabels(app)
    % Compatibility shim for the attached FIONAviewer.m. These labels are
    % bookkeeping only; the visible filename field remains StepsFilename.
    if isempty(app.FileName)
        app.FileName = uilabel(app.LeftPanel, 'Text', '', 'Visible', 'off', ...
            'Position', [1 1 1 1]);
    end
    if isempty(app.FilePath)
        app.FilePath = uilabel(app.LeftPanel, 'Text', '', 'Visible', 'off', ...
            'Position', [1 1 1 1]);
    end
end

function folder = enclosingFolderNamed(startFolder, targetName)
    folder = '';
    candidate = startFolder;
    while ~isempty(candidate)
        [parent, name] = fileparts(candidate);
        if strcmpi(name, targetName)
            folder = candidate;
            return
        end
        if strcmp(parent, candidate), return, end
        candidate = parent;
    end
end
