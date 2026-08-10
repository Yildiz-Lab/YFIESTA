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
    [~, order] = sort({files.name}); files = files(order);
    Report = repmat(struct('InputFile','','OutputFile','','Succeeded',false,'Reason','', ...
        'LongAudit',[],'ShortAudit',[],'Steppiness',[]),numel(files),1);
    RegionScores = cell(numel(files),1);
    for k = 1:numel(files)
        inputFile = fullfile(inputFolder, files(k).name);
        outputFile = fullfile(outputFolder, files(k).name);
        Report(k).InputFile = inputFile; Report(k).OutputFile = outputFile;
        app = [];
        try
            % Existing FIONAviewer loader builds the app and reads data.xy,
            % time, trace, and trace_yx exactly as in interactive use.
            app = FIONAviewer(inputFile);
            if ~opt.ShowViewer, app.fig.Visible = 'off'; end

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
            regions = data.Steppiness.Regions;
            for r = 1:numel(regions), regions(r).OutputFile = outputFile; end
            RegionScores{k} = regions;
            close(app.fig);
        catch ME
            Report(k).Reason = ME.message;
            warning('batchFIONAviewerSIC:TraceFailed', 'Failed for %s: %s', files(k).name, ME.message);
            if ~isempty(app) && isfield(app,'fig') && isvalid(app.fig), close(app.fig); end
        end
    end
    save(fullfile(outputFolder, 'FIONAviewer_SIC5_batch_report.mat'), 'Report', 'opt', '-v7');
    writeFIONASteppinessReports(outputFolder, Report, RegionScores);
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
