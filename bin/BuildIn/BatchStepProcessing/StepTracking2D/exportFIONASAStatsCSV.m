function Report = exportFIONASAStatsCSV(overarchingPath, outputFile)
%EXPORTFIONASASTATSCSV Export S/A FIONAviewer steps and dwell times.
%
%   REPORT = exportFIONASAStatsCSV(PATH) finds the FIONAviewer_SIC output
%   associated with PATH, reads every *_fiona.mat trace in the S and A tier
%   folders, and writes:
%
%       FIONAviewer_SIC/steps_summary_S_A.csv
%
%   PATH may be any of the following:
%       1. The original overarching FIESTA .mat file.
%       2. Its analysis folder.
%       3. The SIC2D_then_FIONAviewerSIC folder.
%       4. The FIONAviewer_SIC folder itself.
%
%   REPORT = exportFIONASAStatsCSV(PATH, OUTPUTFILE) uses the requested CSV
%   filename instead. The CSV layout matches steps_summary_1uM.xlsx: column
%   1 is on-axis step size in nm, column 2 is the preceding dwell time in
%   ms, and one completely blank row separates consecutive molecules.
%
%   Step sizes and dwell boundaries come directly from the final aligned
%   FIONAviewer data.trace array. A changepoint is a row for which column 5
%   is at least one. The dwell paired with that step extends from the prior
%   changepoint (or the first sample) through the sample immediately before
%   the current changepoint.

    if nargin < 1 || isempty(overarchingPath)
        [fileName, folderName] = uigetfile('*.mat', ...
            'Select the overarching FIESTA MAT file');
        if isequal(fileName, 0)
            Report = struct('Cancelled', true);
            return
        end
        overarchingPath = fullfile(folderName, fileName);
    end

    fionaFolder = resolveFIONAFolder(overarchingPath);
    if nargin < 2 || isempty(outputFile)
        outputFile = fullfile(fionaFolder, 'steps_summary_S_A.csv');
    else
        outputFile = char(string(outputFile));
        [outputFolder, ~, extension] = fileparts(outputFile);
        if isempty(outputFolder)
            outputFile = fullfile(fionaFolder, outputFile);
        end
        if isempty(extension)
            outputFile = [outputFile '.csv'];
        elseif ~strcmpi(extension, '.csv')
            error('exportFIONASAStatsCSV:OutputType', ...
                'OUTPUTFILE must have a .csv extension.');
        end
    end

    tierNames = {'S', 'A'};
    files = struct('name', {}, 'folder', {}, 'tier', {}, 'number', {});
    tierCounts = zeros(size(tierNames));
    for t = 1:numel(tierNames)
        tierFolder = fullfile(fionaFolder, tierNames{t});
        if ~isfolder(tierFolder)
            error('exportFIONASAStatsCSV:MissingTier', ...
                'Required tier folder does not exist: %s', tierFolder);
        end
        tierFiles = dir(fullfile(tierFolder, '*_fiona.mat'));
        tierFiles = tierFiles(~startsWith({tierFiles.name}, '._'));
        tierCounts(t) = numel(tierFiles);
        for k = 1:numel(tierFiles)
            entry.name = tierFiles(k).name;
            entry.folder = tierFiles(k).folder;
            entry.tier = tierNames{t};
            entry.number = moleculeNumber(tierFiles(k).name);
            files(end+1) = entry; %#ok<AGROW>
        end
    end

    if isempty(files)
        error('exportFIONASAStatsCSV:NoTraces', ...
            'No S- or A-tier *_fiona.mat traces were found in %s.', fionaFolder);
    end

    [~, order] = sortrows([[files.number].', (1:numel(files)).'], [1 2]);
    files = files(order);

    outputFolder = fileparts(outputFile);
    if ~isempty(outputFolder) && ~isfolder(outputFolder)
        mkdir(outputFolder);
    end
    fileID = fopen(outputFile, 'w');
    if fileID < 0
        error('exportFIONASAStatsCSV:OpenOutput', ...
            'Could not open output file for writing: %s', outputFile);
    end
    cleanup = onCleanup(@() fclose(fileID));
    fprintf(fileID, 'step size (nm),dwell times (ms)\n');

    moleculeNames = strings(numel(files), 1);
    grades = strings(numel(files), 1);
    stepCounts = zeros(numel(files), 1);
    skippedFiles = strings(0, 1);
    writtenMolecules = 0;

    for k = 1:numel(files)
        fullName = fullfile(files(k).folder, files(k).name);
        loaded = load(fullName, 'data');
        if ~isfield(loaded, 'data') || ~isstruct(loaded.data) || ...
                ~isfield(loaded.data, 'trace') || ~isfield(loaded.data, 'time')
            warning('exportFIONASAStatsCSV:InvalidTrace', ...
                'Skipping %s because data.trace or data.time is missing.', fullName);
            skippedFiles(end+1, 1) = string(fullName); %#ok<AGROW>
            continue
        end

        trace = loaded.data.trace;
        time = loaded.data.time(:);
        if ~isnumeric(trace) || size(trace, 2) < 5 || ...
                size(trace, 1) ~= numel(time) || size(trace, 1) < 2
            warning('exportFIONASAStatsCSV:InvalidTraceShape', ...
                'Skipping %s because its trace/time dimensions are invalid.', fullName);
            skippedFiles(end+1, 1) = string(fullName); %#ok<AGROW>
            continue
        end

        changePoints = find(trace(2:end, 5) >= 1) + 1;
        dwellStarts = [1; changePoints(1:end-1)];
        stepSizes = trace(changePoints, 3) - trace(changePoints - 1, 3);
        dwellTimesMs = 1000 .* (time(changePoints - 1) - time(dwellStarts));
        valid = isfinite(stepSizes) & isfinite(dwellTimesMs) & dwellTimesMs >= 0;
        stepSizes = stepSizes(valid);
        dwellTimesMs = dwellTimesMs(valid);

        if writtenMolecules > 0
            fprintf(fileID, '\n');
        end
        for j = 1:numel(stepSizes)
            fprintf(fileID, '%.15g,%.15g\n', stepSizes(j), dwellTimesMs(j));
        end

        writtenMolecules = writtenMolecules + 1;
        moleculeNames(writtenMolecules) = erase(string(files(k).name), '_fiona.mat');
        grades(writtenMolecules) = string(files(k).tier);
        stepCounts(writtenMolecules) = numel(stepSizes);
    end

    moleculeNames = moleculeNames(1:writtenMolecules);
    grades = grades(1:writtenMolecules);
    stepCounts = stepCounts(1:writtenMolecules);

    Report = struct();
    Report.InputPath = char(string(overarchingPath));
    Report.FIONAFolder = fionaFolder;
    Report.OutputFile = outputFile;
    Report.Tiers = tierNames;
    Report.TierFileCounts = tierCounts;
    Report.MoleculeCount = writtenMolecules;
    Report.MoleculeNames = moleculeNames;
    Report.Grades = grades;
    Report.StepCounts = stepCounts;
    Report.TotalStepCount = sum(stepCounts);
    Report.SkippedFiles = skippedFiles;

    fprintf('Exported %d S/A molecules and %d step-dwell pairs to:\n%s\n', ...
        Report.MoleculeCount, Report.TotalStepCount, Report.OutputFile);
end

function fionaFolder = resolveFIONAFolder(inputPath)
    inputPath = char(string(inputPath));
    candidates = strings(0, 1);

    if isfile(inputPath)
        [sourceFolder, baseName] = fileparts(inputPath);
        candidates(end+1) = string(fullfile(sourceFolder, baseName, ...
            'SIC2D_then_FIONAviewerSIC', 'FIONAviewer_SIC'));
        candidates(end+1) = string(fullfile(sourceFolder, ...
            'SIC2D_then_FIONAviewerSIC', 'FIONAviewer_SIC'));
    elseif isfolder(inputPath)
        candidates(end+1) = string(inputPath);
        candidates(end+1) = string(fullfile(inputPath, 'FIONAviewer_SIC'));
        candidates(end+1) = string(fullfile(inputPath, ...
            'SIC2D_then_FIONAviewerSIC', 'FIONAviewer_SIC'));
    else
        error('exportFIONASAStatsCSV:InputNotFound', ...
            'The supplied file or folder does not exist: %s', inputPath);
    end

    for k = 1:numel(candidates)
        candidate = char(candidates(k));
        if isfolder(fullfile(candidate, 'S')) && isfolder(fullfile(candidate, 'A'))
            fionaFolder = candidate;
            return
        end
    end

    error('exportFIONASAStatsCSV:FIONAFolderNotFound', ...
        ['Could not find FIONAviewer_SIC/S and FIONAviewer_SIC/A beneath ' ...
         'the supplied path: %s'], inputPath);
end

function value = moleculeNumber(fileName)
    token = regexp(fileName, '^(-?\d+)', 'tokens', 'once');
    if isempty(token)
        value = Inf;
    else
        value = str2double(token{1});
    end
end
