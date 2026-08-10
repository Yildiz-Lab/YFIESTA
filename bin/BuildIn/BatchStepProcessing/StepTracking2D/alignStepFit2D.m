function [StepStatistics, Alignment] = alignStepFit2D(xy, StepStatistics, Options)
%ALIGNSTEPFIT2D Shift each detected step to the strongest local mean change.
%   In path-constrained mode this is the longitudinal analogue of
%   FIONAviewer's AlignToMaxDx: each decimated SIC boundary is moved to the
%   largest raw-data mean displacement within SearchWindow samples.

    cp = StepStatistics.ChangePoints(:);
    originalCP = cp;
    N = size(xy,1);
    retained = false(size(cp));
    scores = nan(size(cp));
    shifted = nan(size(cp));
    minDwell = Options.MinPointsInDwell;
    [alignmentData, axisAware] = alignmentCoordinate(xy, StepStatistics);
    for j = numel(cp):-1:1
        bounds = [0; cp; N];
        a = bounds(j); b = bounds(j+1); c = bounds(j+2);
        candidates = max(a + minDwell, b - Options.SearchWindow) : ...
                     min(c - minDwell, b + Options.SearchWindow);
        if isempty(candidates), continue, end
        bestScore = -Inf; bestCP = b;
        for candidate = candidates
            left = alignmentData(max(a+1,candidate-Options.SearchWindow+1):candidate,:);
            right = alignmentData(candidate+1:min(c,candidate+Options.SearchWindow),:);
            if axisAware
                score = abs(mean(right) - mean(left));
            else
                score = norm(mean(right,1) - mean(left,1));
            end
            if score > bestScore
                bestScore = score; bestCP = candidate;
            end
        end
        scores(j) = bestScore;
        if bestScore >= Options.MinDisplacement
            retained(j) = true;
            shifted(j) = bestCP;
            cp(j) = bestCP;
        end
    end
    cp = unique(sort(cp(retained)));
    if axisAware
        % Alignment occurs on the unbinned data and can occasionally turn a
        % valid decimated dwell into a brief transient. Remove only a
        % boundary whose removal improves the global 1-D SIC; no dwell-time
        % threshold is imposed here.
        [cp, pruning] = pruneLongitudinalBIC(cp, alignmentData, Options);
    else
        pruning.Method = 'Not applied (joint 2-D segmentation)';
        pruning.RemovedChangePoints = [];
        pruning.FinalSIC = NaN;
    end
    breaks = [0; cp; N];
    [fitXY, centres] = refitDwells(xy, breaks);
    vectors = diff(centres,1,1);
    StepStatistics.ChangePoints = cp;
    StepStatistics.StepFit = fitXY;
    StepStatistics.DwellPositions = centres;
    StepStatistics.NumberOfStepsFound = numel(cp);
    StepStatistics.StepSizeStats = vectors;
    StepStatistics.StepMagnitudes = hypot(vectors(:,1), vectors(:,2));
    StepStatistics.StepAngles = atan2(vectors(:,2), vectors(:,1));
    StepStatistics.StepLengths = diff(breaks);
    StepStatistics.StepsChi2 = localStepReductions(xy, breaks, Options);
    if isfield(StepStatistics,'Options') && isfield(StepStatistics.Options,'DecimationFactor')
        StepStatistics.ChangePointsDecimated = ceil(cp / StepStatistics.Options.DecimationFactor);
        StepStatistics.StepLengthsDecimated = diff([0; ceil(cp / StepStatistics.Options.DecimationFactor); ...
            size(StepStatistics.DecimatedXY,1)]);
    end
    StepStatistics.PostAlignmentSSE = sum(sum((xy-fitXY).^2,2));
    StepStatistics.PostAlignmentSICRaw = rawSIC(xy, breaks, Options);
    if axisAware
        Alignment.Method = 'largest local longitudinal mean displacement + SIC pruning';
    else
        Alignment.Method = 'largest local 2-D mean displacement';
    end
    Alignment.OriginalChangePoints = originalCP;
    Alignment.AlignedChangePoints = cp;
    Alignment.CandidateScores = scores;
    Alignment.ShiftedCandidates = shifted;
    Alignment.BICPruning = pruning;
    Alignment.Options = Options;
    StepStatistics.Alignment = Alignment;
end

function [data, axisAware] = alignmentCoordinate(xy, StepStatistics)
    axisAware = isfield(StepStatistics, 'SegmentationMode') && ...
        strcmpi(char(StepStatistics.SegmentationMode), 'pathConstrained') && ...
        isfield(StepStatistics, 'PathModel') && ~isempty(StepStatistics.PathModel);
    if axisAware
        data = projectStepFit2DPath(xy, StepStatistics.PathModel);
    else
        data = xy;
    end
end

function [cp, Pruning] = pruneLongitudinalBIC(cp, longitudinal, Options)
    % Greedily remove the single boundary with the largest SIC improvement,
    % repeating until every retained boundary has support in the raw
    % longitudinal trace. This is an evidence comparison, not a minimum
    % dwell duration rule.
    N = numel(longitudinal);
    cp = cp(:);
    removed = zeros(0,1);
    [currentSSE, currentSIC] = longitudinalSIC(longitudinal, [0; cp; N], Options);
    while ~isempty(cp)
        bounds = [0; cp; N];
        bestSIC = currentSIC;
        removeIndex = [];
        for j = 1:numel(cp)
            proposed = cp;
            proposed(j) = [];
            [~, proposedSIC] = longitudinalSIC(longitudinal, [0; proposed; N], Options);
            if proposedSIC < bestSIC
                bestSIC = proposedSIC;
                removeIndex = j;
            end
        end
        if isempty(removeIndex), break, end
        removed(end+1,1) = cp(removeIndex); %#ok<AGROW>
        cp(removeIndex) = [];
        [currentSSE, currentSIC] = longitudinalSIC(longitudinal, [0; cp; N], Options);
    end
    Pruning.Method = 'raw longitudinal SIC model comparison';
    Pruning.RemovedChangePoints = removed;
    Pruning.FinalSSE = currentSSE;
    Pruning.FinalSIC = currentSIC;
end

function [totalSSE, value] = longitudinalSIC(data, breaks, Options)
    totalSSE = 0;
    for q = 1:numel(breaks)-1
        rows = breaks(q)+1:breaks(q+1);
        residual = data(rows) - mean(data(rows));
        part = sum(residual.^2);
        if isfield(Options,'CorrectVariance') && Options.CorrectVariance && numel(rows) > 1
            part = part * numel(rows)/(numel(rows)-1);
        end
        totalSSE = totalSSE + part;
    end
    nDwells = numel(breaks)-1;
    N = numel(data);
    value = N * log(max(totalSSE,realmin)/N) + (nDwells+1) * log(N);
end

function reductions = localStepReductions(xy, breaks, Options)
    reductions = zeros(numel(breaks)-2,1);
    for j = 1:numel(reductions)
        reductions(j) = dwellSSE(xy, breaks(j)+1:breaks(j+2), Options) ...
            - dwellSSE(xy, breaks(j)+1:breaks(j+1), Options) ...
            - dwellSSE(xy, breaks(j+1)+1:breaks(j+2), Options);
    end
end

function value = rawSIC(xy, breaks, Options)
    sse = 0;
    for j = 1:numel(breaks)-1
        sse = sse + dwellSSE(xy, breaks(j)+1:breaks(j+1), Options);
    end
    N = size(xy,1); K = numel(breaks)-1;
    value = 2*N*log(max(sse,realmin)/(2*N)) + (2*K+1)*log(N);
end

function sse = dwellSSE(xy, rows, Options)
    n = numel(rows);
    residual = xy(rows,:) - mean(xy(rows,:),1);
    sse = sum(residual(:).^2);
    if isfield(Options,'CorrectVariance') && Options.CorrectVariance && n > 1
        sse = sse * n/(n-1);
    end
end

function [fitXY, centres] = refitDwells(xy, breaks)
    fitXY = zeros(size(xy)); centres = zeros(numel(breaks)-1,2);
    for j = 1:numel(breaks)-1
        rows = breaks(j)+1:breaks(j+1);
        centres(j,:) = mean(xy(rows,:),1);
        fitXY(rows,:) = repmat(centres(j,:),numel(rows),1);
    end
end
