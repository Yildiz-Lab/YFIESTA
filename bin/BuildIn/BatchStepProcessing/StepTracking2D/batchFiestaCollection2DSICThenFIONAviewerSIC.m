function [CollectionReport, CombinedManifest] = batchFiestaCollection2DSICThenFIONAviewerSIC(sourceFiles, outputParent, varargin)
%BATCHFIESTACOLLECTION2DSICTHENFIONAVIEWERSIC Analyze several FIESTA files.
%   Each source MAT file is processed independently with the established
%   2-D SIC -> FIONAviewer SIC pipeline. Final FIONAviewer outputs are also
%   collected into one grade-organized folder using dataset-prefixed names
%   so repeated molecule numbers cannot overwrite one another.

    p = inputParser;
    p.addRequired('sourceFiles', @(x) iscell(x) || isstring(x) || ischar(x));
    p.addRequired('outputParent', @(x) ischar(x) || isstring(x));
    p.addParameter('CopySources', true, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
    p.addParameter('ResumeExisting', true, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
    p.addParameter('MinimumFIONAPreScreenScore', 40, ...
        @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 100);
    p.addParameter('StepFit2DOptions', struct(), @isstruct);
    p.addParameter('FIONAOptions', struct(), @isstruct);
    p.parse(sourceFiles, outputParent, varargin{:});
    opt = p.Results;

    sourceFiles = normalizeFileList(sourceFiles);
    outputParent = char(outputParent);
    if ~isfolder(outputParent), mkdir(outputParent); end
    sourceFolder = fullfile(outputParent, 'Source_FIESTA');
    datasetFolder = fullfile(outputParent, 'Per_Dataset');
    if ~isfolder(sourceFolder), mkdir(sourceFolder); end
    if ~isfolder(datasetFolder), mkdir(datasetFolder); end

    fionaOptions = opt.FIONAOptions;
    if ~isfield(fionaOptions, 'ShowViewer'), fionaOptions.ShowViewer = false; end
    if ~isfield(fionaOptions, 'ResumeExisting'), fionaOptions.ResumeExisting = true; end
    configurationSignature = jsonencode(struct( ...
        'MinimumFIONAPreScreenScore', opt.MinimumFIONAPreScreenScore, ...
        'StepFit2DOptions', opt.StepFit2DOptions, ...
        'FIONAOptions', fionaOptions));

    template = struct('Dataset','','SourceFile','','AnalysisFile','', ...
        'OutputFolder','','Status','Pending','Reason','', ...
        'InputMoleculeCount',0,'PassedMoleculeCount',0, ...
        'Exported2DCount',0,'SelectedFIONACount',0,'FittedFIONACount',0, ...
        'ConfigurationSignature','');
    CollectionReport = repmat(template, numel(sourceFiles), 1);

    for k = 1:numel(sourceFiles)
        sourceFile = sourceFiles{k};
        [~, dataset, extension] = fileparts(sourceFile);
        CollectionReport(k).Dataset = dataset;
        CollectionReport(k).SourceFile = sourceFile;
        CollectionReport(k).OutputFolder = fullfile(datasetFolder, dataset);
        CollectionReport(k).ConfigurationSignature = configurationSignature;
        fprintf('\nCollection dataset %d/%d: %s\n', k, numel(sourceFiles), dataset);

        try
            if ~isfile(sourceFile)
                error('batchFiestaCollection:MissingSource', 'Source file does not exist: %s', sourceFile);
            end
            analysisFile = sourceFile;
            if opt.CopySources
                analysisFile = fullfile(sourceFolder, [dataset extension]);
                copySourceIfNeeded(sourceFile, analysisFile);
            end
            CollectionReport(k).AnalysisFile = analysisFile;

            markerFile = fullfile(CollectionReport(k).OutputFolder, 'collection_dataset_complete.mat');
            if opt.ResumeExisting && isfile(markerFile)
                marker = load(markerFile, 'DatasetCompletion');
                if isfield(marker, 'DatasetCompletion')
                    savedCompletion = conformCompletion(marker.DatasetCompletion, template);
                    if strcmp(savedCompletion.ConfigurationSignature, configurationSignature)
                        CollectionReport(k) = savedCompletion;
                        CollectionReport(k).Status = 'Reused';
                        fprintf('Reusing completed dataset %s.\n', dataset);
                        saveCollectionProgress(outputParent, CollectionReport, opt);
                        continue
                    end
                    fprintf('Analysis settings changed; reprocessing dataset %s.\n', dataset);
                end
            end

            loaded = load(analysisFile, 'Molecule', 'Filament');
            if ~isfield(loaded, 'Molecule') || ~isstruct(loaded.Molecule)
                error('batchFiestaCollection:NoMolecule', ...
                    'File does not contain a struct variable named Molecule: %s', analysisFile);
            end
            CollectionReport(k).InputMoleculeCount = numel(loaded.Molecule);

            args = {'OutputFolder', CollectionReport(k).OutputFolder, ...
                'StepFit2DOptions', opt.StepFit2DOptions, ...
                'FIONAOptions', fionaOptions, ...
                'MinimumFIONAPreScreenScore', opt.MinimumFIONAPreScreenScore};
            if isfield(loaded, 'Filament')
                args = [args {'Filament', loaded.Filament}]; %#ok<AGROW>
            end
            [~, passed, ~, ~, exportReport, fionaReport] = ...
                batchMolecule2DSICThenFIONAviewerSIC(loaded.Molecule, args{:});
            CollectionReport(k).PassedMoleculeCount = numel(passed);
            CollectionReport(k).Exported2DCount = countSucceeded(exportReport);
            CollectionReport(k).SelectedFIONACount = countSelected(CollectionReport(k).OutputFolder);
            CollectionReport(k).FittedFIONACount = countSucceeded(fionaReport);

            if isempty(fionaReport) || all([fionaReport.Succeeded])
                CollectionReport(k).Status = 'Complete';
                CollectionReport(k).Reason = 'Completed without an uncaught trace failure';
                DatasetCompletion = CollectionReport(k); %#ok<NASGU>
                if ~isfolder(CollectionReport(k).OutputFolder)
                    mkdir(CollectionReport(k).OutputFolder);
                end
                save(markerFile, 'DatasetCompletion', '-v7');
            else
                CollectionReport(k).Status = 'Partial';
                CollectionReport(k).Reason = sprintf('%d FIONAviewer trace(s) failed', ...
                    sum(~[fionaReport.Succeeded]));
            end
        catch ME
            CollectionReport(k).Status = 'Failed';
            CollectionReport(k).Reason = ME.message;
            warning('batchFiestaCollection:DatasetFailed', ...
                'Dataset %s failed and the collection will continue: %s', dataset, ME.message);
        end
        saveCollectionProgress(outputParent, CollectionReport, opt);
    end

    CombinedManifest = collectFinalFIONAOutputs(CollectionReport, outputParent);
    save(fullfile(outputParent, 'collection_analysis_report.mat'), ...
        'CollectionReport', 'CombinedManifest', 'opt', '-v7');
    writeCollectionTables(outputParent, CollectionReport, CombinedManifest);
    fprintf('\nCollection complete: %d datasets; %d combined fitted traces.\n', ...
        numel(CollectionReport), numel(CombinedManifest));
end

function files = normalizeFileList(values)
    if ischar(values), files = {values};
    elseif isstring(values), files = cellstr(values(:));
    else, files = values(:);
    end
    files = cellfun(@char, files, 'UniformOutput', false);
end

function copySourceIfNeeded(source, destination)
    shouldCopy = true;
    if isfile(destination)
        sourceInfo = dir(source); destinationInfo = dir(destination);
        shouldCopy = sourceInfo.bytes ~= destinationInfo.bytes;
    end
    if shouldCopy
        fprintf('Copying source file %s...\n', source);
        [ok, message] = copyfile(source, destination, 'f');
        if ~ok
            error('batchFiestaCollection:CopyFailed', ...
                'Could not copy %s to %s: %s', source, destination, message);
        end
    end
end

function completion = conformCompletion(saved, template)
    completion = template;
    names = intersect(fieldnames(saved), fieldnames(template));
    for k = 1:numel(names), completion.(names{k}) = saved.(names{k}); end
end

function count = countSucceeded(report)
    if isempty(report), count = 0;
    elseif isfield(report, 'Succeeded'), count = sum([report.Succeeded]);
    else, count = numel(report);
    end
end

function count = countSelected(outputFolder)
    count = 0;
    reportFile = fullfile(outputFolder, 'combined_SIC_pipeline_report.mat');
    if ~isfile(reportFile), return, end
    values = load(reportFile, 'FIONAPreScreenReport');
    if isfield(values, 'FIONAPreScreenReport') && ...
            ~isempty(values.FIONAPreScreenReport)
        count = sum([values.FIONAPreScreenReport.Selected]);
    end
end

function saveCollectionProgress(outputParent, CollectionReport, opt)
    save(fullfile(outputParent, 'collection_analysis_progress.mat'), ...
        'CollectionReport', 'opt', '-v7');
    writeDatasetTable(outputParent, CollectionReport);
end

function Manifest = collectFinalFIONAOutputs(CollectionReport, outputParent)
    combinedFolder = fullfile(outputParent, 'Combined_FIONAviewer_SIC');
    tiers = {'S','A','B','C_or_lower','all'};
    for k = 1:numel(tiers)
        folder = fullfile(combinedFolder, tiers{k});
        if ~isfolder(folder), mkdir(folder); end
    end
    template = struct('Dataset','','OriginalName','','Grade','','Score',NaN, ...
        'SourceMAT','','CombinedMAT','','CombinedFigure','','CombinedFlagFigure','');
    Manifest = repmat(template, 0, 1);

    for d = 1:numel(CollectionReport)
        sourceRoot = fullfile(CollectionReport(d).OutputFolder, 'FIONAviewer_SIC');
        if ~isfolder(sourceRoot), continue, end
        for t = 1:4
            sourceTier = fullfile(sourceRoot, tiers{t});
            files = dir(fullfile(sourceTier, '*_fiona.mat'));
            files = files(~startsWith({files.name}, '._'));
            for k = 1:numel(files)
                sourceMAT = fullfile(files(k).folder, files(k).name);
                loaded = load(sourceMAT, 'data');
                if ~isfield(loaded, 'data') || ~isfield(loaded.data, 'Steppiness')
                    warning('batchFiestaCollection:UnscoredTrace', ...
                        'Skipping combined copy of unscored trace %s.', sourceMAT);
                    continue
                end
                grade = upper(char(loaded.data.Steppiness.Trace.Grade));
                gradeFolder = finalGradeFolder(grade);
                [~, originalName] = fileparts(sourceMAT);
                combinedBase = [safeFilename(CollectionReport(d).Dataset) '__' originalName];
                destinationFolder = fullfile(combinedFolder, gradeFolder);
                destinationMAT = fullfile(destinationFolder, [combinedBase '.mat']);
                checkedCopy(sourceMAT, destinationMAT);
                checkedCopy(sourceMAT, fullfile(combinedFolder, 'all', [combinedBase '.mat']));

                sourceFigure = fullfile(files(k).folder, [originalName '.fig']);
                destinationFigure = fullfile(destinationFolder, [combinedBase '.fig']);
                if isfile(sourceFigure), checkedCopy(sourceFigure, destinationFigure); end
                sourceFlag = fullfile(files(k).folder, [originalName '_review_flags.fig']);
                destinationFlag = fullfile(destinationFolder, [combinedBase '_review_flags.fig']);
                if isfile(sourceFlag), checkedCopy(sourceFlag, destinationFlag); end

                entry = template;
                entry.Dataset = CollectionReport(d).Dataset;
                entry.OriginalName = originalName;
                entry.Grade = grade;
                entry.Score = loaded.data.Steppiness.Trace.Score;
                entry.SourceMAT = sourceMAT;
                entry.CombinedMAT = destinationMAT;
                if isfile(destinationFigure), entry.CombinedFigure = destinationFigure; end
                if isfile(destinationFlag), entry.CombinedFlagFigure = destinationFlag; end
                Manifest(end+1,1) = entry; %#ok<AGROW>
            end
        end
    end
end

function checkedCopy(source, destination)
    [ok, message] = copyfile(source, destination, 'f');
    if ~ok
        error('batchFiestaCollection:CombinedCopyFailed', ...
            'Could not copy %s to %s: %s', source, destination, message);
    end
end

function value = safeFilename(value)
    value = regexprep(char(value), '[^A-Za-z0-9._-]', '_');
end

function folder = finalGradeFolder(grade)
    switch upper(char(grade))
        case 'S', folder = 'S';
        case 'A', folder = 'A';
        case 'B', folder = 'B';
        otherwise, folder = 'C_or_lower';
    end
end

function writeCollectionTables(outputParent, CollectionReport, Manifest)
    writeDatasetTable(outputParent, CollectionReport);
    if isempty(Manifest)
        writetable(struct2table(repmat(struct('Dataset','','OriginalName','', ...
            'Grade','','Score',NaN,'SourceMAT','','CombinedMAT','', ...
            'CombinedFigure','','CombinedFlagFigure',''), 0, 1)), ...
            fullfile(outputParent, 'combined_trace_manifest.csv'));
    else
        writetable(struct2table(Manifest), ...
            fullfile(outputParent, 'combined_trace_manifest.csv'));
    end
end

function writeDatasetTable(outputParent, CollectionReport)
    writetable(struct2table(CollectionReport), ...
        fullfile(outputParent, 'dataset_processing_report.csv'));
end
