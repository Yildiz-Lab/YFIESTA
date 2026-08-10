function [AllMolecules, PassedMolecules, StepFit2DReport, Options, FIONAExportReport, FIONAReport] = ...
    batchMolecule2DSICThenFIONAviewerSIC(Molecule, varargin)
%BATCHMOLECULE2DSICTHENFIONAVIEWERSIC Complete 2-D SIC -> FIONA SIC workflow.
%   1) uses processMoleculeStepFit2D for QC, EFO truncation/flags, 2-D SIC,
%      and alignment;
%   2) exports passed traces to a FIONA_Input_2DSIC subfolder;
%   3) invokes the existing folder-based batchFIONAviewerSIC on that folder;
%   4) saves FIONAviewer SIC traces, figures, and revised review reports.

    thisFolder = fileparts(mfilename('fullpath'));
    fionaBatchFolder = fullfile(thisFolder, 'FIONAviewerBatchSIC');
    % The standalone package keeps FIONAviewerBatchSIC beside this file,
    % whereas YFIESTA places both folders directly under BuildIn.
    if ~isfolder(fionaBatchFolder)
        fionaBatchFolder = fullfile(fileparts(thisFolder), 'FIONAviewerBatchSIC');
    end
    if ~isfolder(fionaBatchFolder)
        error('batchMolecule2DSICThenFIONAviewerSIC:FIONABatchPath', ...
            'Could not locate the companion FIONAviewerBatchSIC folder.');
    end
    addpath(fionaBatchFolder, '-begin');
    rehash;
    expectedBatchFunction = fullfile(fionaBatchFolder, 'batchFIONAviewerSIC.m');
    resolvedBatchFunction = which('batchFIONAviewerSIC');
    if ~strcmpi(resolvedBatchFunction, expectedBatchFunction)
        error('batchMolecule2DSICThenFIONAviewerSIC:FIONABatchResolution', ...
            ['MATLAB resolved batchFIONAviewerSIC to %s instead of the intended ' ...
             'implementation at %s. Run clear functions, then retry.'], ...
            resolvedBatchFunction, expectedBatchFunction);
    end
    p = inputParser;
    p.addRequired('Molecule', @isstruct);
    p.addParameter('OutputFolder', '', @(x) ischar(x) || isstring(x));
    p.addParameter('OutputFolderName', 'SIC2D_then_FIONAviewerSIC', @(x) ischar(x) || isstring(x));
    p.addParameter('StepFit2DOptions', struct(), @isstruct);
    p.addParameter('FIONAOptions', struct(), @isstruct);
    p.addParameter('Filament', [], @(x) isempty(x) || isstruct(x));
    % Retain S through D by default. F is still exported from 2-D SIC but
    % does not enter the slower FIONAviewer SIC stage.
    p.addParameter('MinimumFIONAPreScreenScore', 40, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 100);
    p.parse(Molecule, varargin{:}); opt = p.Results;

    if isempty(opt.OutputFolder)
        outputParent = uigetdir(pwd, 'Select parent folder for combined SIC pipeline');
        if isequal(outputParent,0)
            AllMolecules = Molecule; PassedMolecules = Molecule([]);
            StepFit2DReport = struct([]); Options = stepFit2DDefaultOptions();
            FIONAExportReport = struct([]); FIONAReport = struct([]); return
        end
        outputFolder = fullfile(outputParent, char(opt.OutputFolderName));
    else
        outputFolder = char(opt.OutputFolder);
    end
    if ~isfolder(outputFolder), mkdir(outputFolder); end

    Options = mergeStepFit2DOptions(stepFit2DDefaultOptions(), opt.StepFit2DOptions);
    [AllMolecules, PassedMolecules, StepFit2DReport, Options] = ...
        processMoleculeStepFit2D(Molecule, Options);
    save(fullfile(outputFolder, 'Molecule_2DSIC_all.mat'), ...
        'AllMolecules', 'StepFit2DReport', 'Options', '-v7.3');
    save(fullfile(outputFolder, 'Molecule_2DSIC_passed_truncated.mat'), ...
        'PassedMolecules', 'StepFit2DReport', 'Options', '-v7.3');
    % This separate file uses the exact Molecule/Filament variable names,
    % schema, field order, and row-vector orientation required by FIESTA's
    % LoadTracks. It contains only the passed, EFO-truncated molecules.
    saveStepFit2DForFiesta(fullfile(outputFolder, 'Molecule_2DSIC_passed_FIESTA.mat'), ...
        PassedMolecules, opt.Filament);

    fionaInputFolder = fullfile(outputFolder, 'FIONA_Input_2DSIC');
    FIONAExportReport = exportStepFit2DToFIONA(PassedMolecules, fionaInputFolder);
    try
        TwoDFigureReport = export2DSICFigures(PassedMolecules, FIONAExportReport, fionaInputFolder);
    catch figureME
        warning('batchMolecule2DSICThenFIONAviewerSIC:FigureExport', ...
            '2-D figure export failed, but FIONAviewer SIC will continue: %s', figureME.message);
        TwoDFigureReport = emptyTwoDFigureReport(numel(PassedMolecules), figureME.message);
    end
    [FIONAPreScreenReport, selectedInputFolder, FIONAExportReport, TwoDFigureReport] = preScreenFIONAInputs( ...
        FIONAExportReport, TwoDFigureReport, fionaInputFolder, Options.Fit.DecimationFactor, ...
        opt.MinimumFIONAPreScreenScore);
    save(fullfile(outputFolder, 'FIONA_Input_2DSIC_export_report.mat'), ...
        'FIONAExportReport', 'TwoDFigureReport', 'FIONAPreScreenReport', '-v7');
    exported = FIONAExportReport([FIONAExportReport.Succeeded]);
    if isempty(exported)
        FIONAReport = struct([]);
        warning('batchMolecule2DSICThenFIONAviewerSIC:NoExports', ...
            'No passed 2-D SIC molecules were exported; FIONAviewer SIC was not run.');
        return
    end
    selected = FIONAPreScreenReport([FIONAPreScreenReport.Selected]);
    fprintf('FIONAviewer pre-screen: %d of %d exported traces met the score threshold (score >= %.0f).\n', ...
        numel(selected), numel(exported), opt.MinimumFIONAPreScreenScore);
    if isempty(selected)
        FIONAReport = struct([]);
        save(fullfile(outputFolder, 'combined_SIC_pipeline_report.mat'), ...
            'StepFit2DReport', 'FIONAExportReport', 'TwoDFigureReport', ...
            'FIONAPreScreenReport', 'FIONAReport', 'Options', 'opt', '-v7');
        warning('batchMolecule2DSICThenFIONAviewerSIC:NoSelectedTraces', ...
            'No traces met the FIONAviewer pre-screen criterion; FIONAviewer SIC was not run.');
        return
    end

    fionaOutputFolder = fullfile(outputFolder, 'FIONAviewer_SIC');
    fionaArgs = structToNameValue(opt.FIONAOptions);
    FIONAReport = batchFIONAviewerSIC(selectedInputFolder, ...
        fionaArgs{:}, 'OutputFolder', fionaOutputFolder);
    FIONACompletion = summarizeFIONAOutputs(selected, FIONAReport);
    save(fullfile(outputFolder, 'combined_SIC_pipeline_report.mat'), ...
        'StepFit2DReport', 'FIONAExportReport', 'TwoDFigureReport', ...
        'FIONAPreScreenReport', 'FIONAReport', 'FIONACompletion', ...
        'Options', 'opt', '-v7');
    if FIONACompletion.FittedTraceCount == 0
        error('batchMolecule2DSICThenFIONAviewerSIC:NoFIONAOutputs', ...
            ['The 2-D SIC traces were exported, but FIONAviewer produced no final traces. ' ...
             'See FIONAviewer_SIC/FIONAviewer_SIC5_batch_report.mat for the first failure.']);
    end
    if FIONACompletion.FittedTraceCount < FIONACompletion.ExportedTraceCount
        warning('batchMolecule2DSICThenFIONAviewerSIC:PartialFIONAOutput', ...
            'Created %d of %d final FIONA traces. See the combined report for failures.', ...
            FIONACompletion.FittedTraceCount, FIONACompletion.ExportedTraceCount);
    end
    fprintf(['FIONA exports: %d input traces; %d 2-D figure pairs; %d fitted traces; ' ...
        '%d standard figures; %d review-flag figures.\n'], FIONACompletion.ExportedTraceCount, ...
        sum([TwoDFigureReport.Succeeded]), ...
        FIONACompletion.FittedTraceCount, FIONACompletion.FigureCount, ...
        FIONACompletion.FlagFigureCount);
end

function args = structToNameValue(values)
    names = fieldnames(values); args = cell(1,2*numel(names));
    for k = 1:numel(names)
        args{2*k-1} = names{k}; args{2*k} = values.(names{k});
    end
end

function Summary = summarizeFIONAOutputs(ExportReport, FIONAReport)
    if isfield(ExportReport, 'Succeeded')
        Summary.ExportedTraceCount = sum([ExportReport.Succeeded]);
    elseif isfield(ExportReport, 'Selected')
        Summary.ExportedTraceCount = sum([ExportReport.Selected]);
    else
        Summary.ExportedTraceCount = numel(ExportReport);
    end
    successful = FIONAReport([FIONAReport.Succeeded]);
    Summary.FittedTraceCount = numel(successful);
    Summary.OutputTraceFiles = {successful.OutputFile};
    if isempty(successful)
        Summary.FigureFiles = {};
        Summary.FlagFigureFiles = {};
        Summary.FigureCount = 0;
        Summary.FlagFigureCount = 0;
        return
    end
    figureFiles = cell(size(successful));
    flagFigureFiles = cell(size(successful));
    for k = 1:numel(successful)
        [folder, name] = fileparts(successful(k).OutputFile);
        figureFiles{k} = fullfile(folder, [name '.fig']);
        if isfield(successful, 'FlagFigure') && ~isempty(successful(k).FlagFigure)
            flagFigureFiles{k} = successful(k).FlagFigure;
        else
            % Support an older batchFIONAviewerSIC report while still
            % recognizing a review figure written beside the trace.
            flagFigureFiles{k} = fullfile(folder, [name '_review_flags.fig']);
        end
    end
    Summary.FigureFiles = figureFiles;
    Summary.FlagFigureFiles = flagFigureFiles;
    Summary.FigureCount = sum(cellfun(@isfile, figureFiles));
    Summary.FlagFigureCount = sum(cellfun(@isfile, flagFigureFiles));
end

function Report = export2DSICFigures(Molecules, ExportReport, outputFolder)
    Report = repmat(struct('Name','','TrajectoryFile','','ProjectionFile', ...
        '','Succeeded',false,'Reason',''), numel(Molecules), 1);
    for k = 1:numel(Molecules)
        Report(k).Name = ExportReport(k).Name;
        if ~ExportReport(k).Succeeded
            Report(k).Reason = 'FIONA trace export failed';
            continue
        end
        trajectoryFigure = [];
        projectionFigure = [];
        try
            fprintf('Exporting figure %s (2-D SIC trajectory and projection)...\n', Report(k).Name);
            results = Molecules(k).Results;
            Figures = plotSICstepFinder2D(results(:,3:4), Molecules(k).StepFit2D, ...
                'Time', results(:,2), 'PathModel', Molecules(k).PathFit2D, ...
                'MoleculeName', figureMoleculeName(Molecules(k), k), 'Visible', 'off');
            trajectoryFigure = Figures.Trajectory;
            projectionFigure = Figures.Projection;
            [~, stem] = fileparts(ExportReport(k).File);
            Report(k).TrajectoryFile = fullfile(outputFolder, [stem '_2DSIC_trajectory.fig']);
            Report(k).ProjectionFile = fullfile(outputFolder, [stem '_2DSIC_projection.fig']);
            % These figures are created hidden so the batch can run without
            % opening many windows.  Persist them as visible, however:
            % otherwise opening the saved FIG loads it successfully but
            % shows no window.  The temporary batch figures are hidden and
            % closed again immediately after writing.
            trajectoryFigure.Visible = 'on';
            projectionFigure.Visible = 'on';
            % Match FIONAviewer's portable FIG export path.  drawnow
            % resolves linked axes, legends, and xline objects before the
            % figure is serialized; saveas writes the conventional FIG
            % representation used by FIONAviewer's ExportFigure helper.
            drawnow;
            saveas(trajectoryFigure, Report(k).TrajectoryFile, 'fig');
            saveas(projectionFigure, Report(k).ProjectionFile, 'fig');
            trajectoryFigure.Visible = 'off';
            projectionFigure.Visible = 'off';
            Report(k).Succeeded = true;
            Report(k).Reason = 'Exported';
        catch ME
            Report(k).Reason = ME.message;
            warning('batchMolecule2DSICThenFIONAviewerSIC:FigureExport', ...
                'Could not export 2-D SIC figures for %s: %s', Report(k).Name, ME.message);
        end
        if ~isempty(trajectoryFigure) && isvalid(trajectoryFigure), close(trajectoryFigure); end
        if ~isempty(projectionFigure) && isvalid(projectionFigure), close(projectionFigure); end
    end
end

function [Report, selectedFolder, ExportReport, FigureReport] = preScreenFIONAInputs(ExportReport, FigureReport, inputFolder, downsampleFactor, minimumScore)
    % Category folders retain every 2-D SIC FIONA trace for review. The
    % separate staging folder lets the existing non-recursive FIONAviewer
    % batch process only grades that meet the requested score threshold.
    selectedFolder = fullfile(inputFolder, 'For_FIONAviewer_SIC');
    if ~isfolder(selectedFolder), mkdir(selectedFolder); end
    Report = repmat(struct('Name','','SourceFile','','Score',NaN,'Grade','', ...
        'GradeFolder','','CategorizedFile','','TrajectoryFigure','', ...
        'ProjectionFigure','','SelectedFile','', ...
        'Selected',false,'Reason',''), numel(ExportReport), 1);
    for k = 1:numel(ExportReport)
        Report(k).Name = ExportReport(k).Name;
        Report(k).SourceFile = ExportReport(k).File;
        if ~ExportReport(k).Succeeded
            Report(k).Reason = '2-D FIONA export failed';
            continue
        end
        try
            loaded = load(ExportReport(k).File, 'data');
            data = loaded.data;
            data.Steppiness2DPreScreen = scoreFIONASteppiness(data, ...
                'DownsampleFactor', downsampleFactor);
            Report(k).Score = data.Steppiness2DPreScreen.Trace.Score;
            Report(k).Grade = data.Steppiness2DPreScreen.Trace.Grade;
            Report(k).Selected = Report(k).Score >= minimumScore;
            gradeFolder = fullfile(inputFolder, gradeFolderName(Report(k).Grade));
            if ~isfolder(gradeFolder), mkdir(gradeFolder); end
            [~, filename, extension] = fileparts(ExportReport(k).File);
            categorizedFile = fullfile(gradeFolder, [filename extension]);
            save(categorizedFile, 'data', '-v7');
            % The root export is replaced by its categorized copy, so all
            % per-trace MAT files are organized only by their 2-D grade.
            if ~strcmp(ExportReport(k).File, categorizedFile) && isfile(ExportReport(k).File)
                delete(ExportReport(k).File);
            end
            ExportReport(k).File = categorizedFile;
            ExportReport(k).Reason = ['Exported to ' gradeFolderName(Report(k).Grade)];
            Report(k).GradeFolder = gradeFolder;
            Report(k).CategorizedFile = categorizedFile;
            [FigureReport, Report(k)] = move2DSICFiguresToGrade( ...
                FigureReport, k, gradeFolder, Report(k));
            if Report(k).Selected
                Report(k).SelectedFile = fullfile(selectedFolder, [filename extension]);
                save(Report(k).SelectedFile, 'data', '-v7');
                Report(k).Reason = 'Selected for FIONAviewer SIC';
            else
                Report(k).Reason = 'Below FIONAviewer pre-screen threshold';
            end
        catch ME
            Report(k).Reason = ME.message;
            warning('batchMolecule2DSICThenFIONAviewerSIC:PreScreen', ...
                'Could not score 2-D SIC export %s: %s', Report(k).Name, ME.message);
        end
    end
end

function [FigureReport, preScreenEntry] = move2DSICFiguresToGrade(FigureReport, index, gradeFolder, preScreenEntry)
    % Keep each 2-D SIC FIG pair with its categorized FIONA MAT-file.  The
    % FIONAviewer staging copy contains MAT files only, so this relocation
    % does not affect the downstream non-recursive batch.
    if index > numel(FigureReport) || ~FigureReport(index).Succeeded
        return
    end
    [FigureReport(index).TrajectoryFile, preScreenEntry.TrajectoryFigure] = ...
        moveFigureToFolder(FigureReport(index).TrajectoryFile, gradeFolder);
    [FigureReport(index).ProjectionFile, preScreenEntry.ProjectionFigure] = ...
        moveFigureToFolder(FigureReport(index).ProjectionFile, gradeFolder);
end

function [destination, recordedDestination] = moveFigureToFolder(source, destinationFolder)
    destination = source;
    recordedDestination = '';
    if isempty(source) || ~isfile(source), return, end
    [~, name, extension] = fileparts(source);
    destination = fullfile(destinationFolder, [name extension]);
    if ~strcmp(source, destination)
        [moved, message] = movefile(source, destination, 'f');
        if ~moved
            warning('batchMolecule2DSICThenFIONAviewerSIC:FigureOrganize', ...
                'Could not move %s into %s: %s', source, destinationFolder, message);
            destination = source;
        end
    end
    recordedDestination = destination;
end

function folder = gradeFolderName(grade)
    switch upper(char(grade))
        case 'S', folder = 'S';
        case 'A', folder = 'A';
        case 'B', folder = 'B';
        otherwise, folder = 'C_or_lower';
    end
end

function Report = emptyTwoDFigureReport(count, reason)
    Report = repmat(struct('Name','','TrajectoryFile','','ProjectionFile','', ...
        'Succeeded',false,'Reason',reason), count, 1);
end

function name = figureMoleculeName(Molecule, index)
    name = sprintf('Molecule %d', index);
    if ~isfield(Molecule, 'Name') || isempty(Molecule.Name), return, end
    value = Molecule.Name;
    if iscell(value), value = value{1}; end
    if isnumeric(value), value = num2str(value); elseif isstring(value), value = char(value); end
    value = strtrim(value);
    if isempty(value), return, end
    if strncmpi(value, 'molecule', numel('molecule'))
        name = value;
    else
        name = ['Molecule ' value];
    end
end
