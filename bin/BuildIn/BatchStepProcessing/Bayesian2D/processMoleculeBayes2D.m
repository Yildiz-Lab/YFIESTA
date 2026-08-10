function [AllMolecules, PassedMolecules, Report, Options] = processMoleculeBayes2D(Molecule, UserOptions)
%PROCESSMOLECULEBAYES2D Raw Molecule to globally segmented Bayesian 2-D fits.
    if nargin < 2, UserOptions = struct; end
    Options = mergeStepFit2DOptions(bayesStepFit2DDefaultOptions(), UserOptions);
    AllMolecules = Molecule;
    for i = 1:numel(AllMolecules)
        AllMolecules(i).StepFit2D = [];
        AllMolecules(i).PathFit2D = [];
        AllMolecules(i).StepFit2DStatus = [];
    end
    PassedMolecules = AllMolecules([]);
    Report = repmat(struct('Index',[],'Name','','Passed',false,'Reason','','Quality',[],'EFO',[]),numel(Molecule),1);

    for i = 1:numel(Molecule)
        label = moleculeLabel(Molecule(i),i);
        Report(i).Index = i; Report(i).Name = label;
        status.Stage = 'quality filter'; status.Passed = false; status.Reason = '';
        [qualityPassed,quality] = assessStepFit2DMolecule(Molecule(i),Options.Filter);
        Report(i).Quality = quality;
        if Options.Filter.Enabled && ~qualityPassed
            status.Reason = quality.Reason; AllMolecules(i).StepFit2DStatus = status; Report(i).Reason = status.Reason; continue
        end
        candidate = AllMolecules(i);
        [candidate,efo] = truncateHighEFOStepFit2D(candidate,Options.EFO);
        Report(i).EFO = efo;
        if ~efo.Passed
            status.Stage = 'EFO truncation'; status.Reason = efo.Reason; AllMolecules(i).StepFit2DStatus = status; Report(i).Reason = status.Reason; continue
        end
        if size(candidate.Results,1) < Options.Filter.MinPoints
            status.Stage = 'post-EFO quality check'; status.Reason = 'Too short after EFO truncation';
            AllMolecules(i).StepFit2DStatus = status; Report(i).Reason = status.Reason; continue
        end
        xy = double(candidate.Results(:,3:4)); time = double(candidate.Results(:,2));
        if Options.Plot.EFOAfterTruncation
            plotStepFit2DEFO(time,candidate.Results(:,Options.EFO.Column),efo,label);
        end
        pathModel = fitStepFit2DPath(xy,Options.Path.Degree);
        covariance = [];
        if Options.Bayes.UsePathWeightedNoise
            covariance = pathWeightedNoiseCovariance(xy,pathModel);
        end
        fit = BayesStepFinder2D(xy,'MinPointsInStep',Options.Bayes.MinPointsInStep, ...
            'DecimationFactor',Options.Bayes.DecimationFactor, ...
            'ExpectedDwellPoints',Options.Bayes.ExpectedDwellPoints, ...
            'PriorMeanStrength',Options.Bayes.PriorMeanStrength, ...
            'NoiseCovariance',covariance);
        if Options.Alignment.Enabled
            [fit,alignment] = alignStepFit2D(xy,fit,Options.Alignment); fit.Alignment = alignment;
        end
        fit.PathModel = pathModel;
        fit.EFOTruncation = efo;
        if Options.Plot.FinalFits
            plotSICstepFinder2D(xy,fit,'Time',time,'PathModel',pathModel, ...
                'MoleculeName',label,'Units',Options.Plot.Units);
        end
        status.Stage = 'complete'; status.Passed = true; status.Reason = 'Passed';
        candidate.StepFit2D = fit; candidate.PathFit2D = pathModel; candidate.StepFit2DStatus = status;
        AllMolecules(i) = candidate; PassedMolecules(end+1) = candidate; %#ok<AGROW>
        Report(i).Passed = true; Report(i).Reason = 'Passed';
    end
end

function covariance = pathWeightedNoiseCovariance(xy,pathModel)
    [on,off] = projectStepFit2DPath(xy,pathModel);
    scales = 1.4826 * median(abs(diff([on off])-median(diff([on off]),1)),1) / sqrt(2);
    scales = max(scales,eps(max(abs(xy(:)))));
    u = pathModel.AxisU(:); v = pathModel.AxisV(:);
    covariance = scales(1)^2*(u*u.') + scales(2)^2*(v*v.');
end

function label = moleculeLabel(Molecule,index)
    label = sprintf('Molecule %d',index);
    if isfield(Molecule,'Name') && ~isempty(Molecule.Name)
        raw = Molecule.Name; if isnumeric(raw), raw = num2str(raw); elseif isstring(raw), raw = char(raw); end
        raw = strtrim(raw); if ~isempty(raw)
            if strncmpi(raw,'molecule',8), label = raw; else, label = ['Molecule ' raw]; end
        end
    end
end
