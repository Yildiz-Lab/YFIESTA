function [AllMolecules, PassedMolecules, Report, Options] = processMoleculeStepFit2D(Molecule, UserOptions)
%PROCESSMOLECULESTEPFIT2D Complete modular raw-Molecule to 2-D-stepfit flow.
%   The output AllMolecules preserves every input molecule and adds:
%       StepFit2DStatus  - pass/fail status and QC/EFO information
%       StepFit2D        - final aligned 2-D fit for passed molecules
%       PathFit2D        - linear or quadratic path model for passed molecules
%   PassedMolecules contains only passed, EFO-truncated molecules.

    if nargin < 2, UserOptions = struct; end
    Options = mergeStepFit2DOptions(stepFit2DDefaultOptions(), UserOptions);
    if ~isstruct(Molecule)
        error('processMoleculeStepFit2D:Input', 'Molecule must be a structure array.');
    end
    AllMolecules = Molecule;
    % Predeclare pipeline fields on every element so passed and rejected
    % molecules remain one compatible MATLAB structure array.
    for i = 1:numel(AllMolecules)
        AllMolecules(i).StepFit2D = [];
        AllMolecules(i).PathFit2D = [];
        AllMolecules(i).StepFit2DStatus = [];
    end
    PassedMolecules = AllMolecules([]);
    Report = repmat(struct('Index',[],'Name','','Passed',false,'Reason','','Quality',[],'EFO',[]), numel(Molecule), 1);

    for i = 1:numel(Molecule)
        Report(i).Index = i;
        moleculeName = makeMoleculeName(Molecule(i), i);
        Report(i).Name = moleculeName;
        status.Stage = 'quality filter'; status.Passed = false; status.Reason = '';
        [qualityPassed, quality] = assessStepFit2DMolecule(Molecule(i), Options.Filter);
        Report(i).Quality = quality;
        if Options.Filter.Enabled && ~qualityPassed
            status.Reason = quality.Reason;
            AllMolecules(i).StepFit2DStatus = status;
            Report(i).Reason = status.Reason;
            continue
        end

        candidate = AllMolecules(i);
        [candidate, efo] = truncateHighEFOStepFit2D(candidate, Options.EFO);
        Report(i).EFO = efo;
        if ~efo.Passed
            status.Stage = 'EFO truncation'; status.Reason = efo.Reason;
            AllMolecules(i).StepFit2DStatus = status;
            Report(i).Reason = status.Reason;
            continue
        end
        if size(candidate.Results,1) < Options.Filter.MinPoints
            status.Stage = 'post-EFO quality check';
            status.Reason = 'Too short after EFO truncation';
            AllMolecules(i).StepFit2DStatus = status;
            Report(i).Reason = status.Reason;
            continue
        end

        xy = candidate.Results(:,3:4);
        time = candidate.Results(:,2);
        if Options.Plot.EFOAfterTruncation
            plotStepFit2DEFO(time, candidate.Results(:,Options.EFO.Column), efo, moleculeName);
        end
        pathModel = fitStepFit2DPath(xy, Options.Path.Degree);
        initialFit = SICstepFinder2D(xy, ...
            'MinPointsInStep', Options.Fit.MinPointsInStep, ...
            'CorrectVariance', Options.Fit.CorrectVariance, ...
            'DecimationFactor', Options.Fit.DecimationFactor, ...
            'MaxSteps', Options.Fit.MaxSteps, ...
            'SegmentationMode', Options.Fit.SegmentationMode, ...
            'PathModel', pathModel);
        % The alignment stage uses this model to seek the strongest local
        % longitudinal change rather than a transient transverse excursion.
        initialFit.PathModel = pathModel;
        if Options.Plot.InitialFits
            plotSICstepFinder2D(xy, initialFit, 'Time', time, 'PathModel', pathModel, ...
                'MoleculeName', moleculeName, 'Units', Options.Plot.Units);
        end

        if Options.Alignment.Enabled
            alignmentOptions = Options.Alignment;
            alignmentOptions.CorrectVariance = Options.Fit.CorrectVariance;
            if isempty(alignmentOptions.SearchWindow)
                % SIC finds a boundary on decimated bins. Refine it only
                % within one such bin-width on either side in raw samples,
                % matching the scale used by FIONAviewer's AlignToMaxDx.
                alignmentOptions.SearchWindow = Options.Fit.DecimationFactor;
                alignmentOptions.SearchWindowSource = 'DecimationFactor';
            else
                alignmentOptions.SearchWindowSource = 'User supplied';
            end
            [finalFit, alignment] = alignStepFit2D(xy, initialFit, alignmentOptions);
        else
            finalFit = initialFit;
            alignment.Method = 'Not applied'; alignment.Options = Options.Alignment;
            finalFit.Alignment = alignment;
        end
        finalFit.PreAlignmentChangePoints = initialFit.ChangePoints;
        finalFit.PreAlignmentDwellPositions = initialFit.DwellPositions;
        finalFit.PathModel = pathModel;
        finalFit.EFOTruncation = efo;
        finalFit.EFOFlags = flagStepFit2DEFORegions(time, ...
            candidate.Results(:,Options.EFO.Column), efo.Threshold, ...
            Options.EFO.FlagMultiplier);
        if Options.Plot.AlignedFits
            plotSICstepFinder2D(xy, finalFit, 'Time', time, 'PathModel', pathModel, ...
                'MoleculeName', moleculeName, 'Units', Options.Plot.Units);
        end

        status.Stage = 'complete'; status.Passed = true; status.Reason = 'Passed';
        candidate.StepFit2D = finalFit;
        candidate.PathFit2D = pathModel;
        candidate.StepFit2DStatus = status;
        AllMolecules(i) = candidate;
        PassedMolecules(end+1) = candidate; %#ok<AGROW>
        Report(i).Passed = true;
        Report(i).Reason = 'Passed';
    end
end

function name = makeMoleculeName(Molecule, index)
    name = sprintf('Molecule %d', index);
    if ~isfield(Molecule, 'Name') || isempty(Molecule.Name)
        return
    end
    if iscell(Molecule.Name)
        rawName = Molecule.Name{1};
    else
        rawName = Molecule.Name;
    end
    if isnumeric(rawName)
        rawName = num2str(rawName);
    elseif isstring(rawName)
        rawName = char(rawName);
    end
    rawName = strtrim(rawName);
    if isempty(rawName), return, end
    if strncmpi(rawName, 'molecule', numel('molecule'))
        name = rawName;
    else
        name = ['Molecule ' rawName];
    end
end
