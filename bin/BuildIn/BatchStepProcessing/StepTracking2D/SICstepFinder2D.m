function StepStatistics = SICstepFinder2D(xy, varargin)
%SICSTEPFINDER2D Find piecewise-constant steps in a 2-D trajectory using SIC.
%
%   S = SICstepFinder2D(XY) segments the N-by-2 coordinates XY into dwells
%   with constant 2-D positions. In joint2D mode it is a greedy 2-D SIC
%   analogue of SICstepFinder. In pathConstrained mode it uses global
%   dynamic-programming optimisation of SIC on longitudinal motion.
%
%   S = SICstepFinder2D(XY, 'MinPointsInStep', M, ...)
%   accepts these name/value options:
%       MinPointsInStep - minimum number of localizations per dwell (2)
%       CorrectVariance - apply n/(n-1) small-sample correction (true)
%       DecimationFactor - average this many consecutive points before
%                          fitting; the returned StepFit is expanded back
%                          to the original point count (1)
%       MaxSteps        - maximum number of accepted steps (Inf)
%       SegmentationMode - 'joint2D' (default) or 'pathConstrained'.
%                          The latter selects changepoints from the
%                          longitudinal coordinate of PathModel only,
%                          then returns full 2-D dwell means.
%       PathModel       - model from fitStepFit2DPath; required for
%                          SegmentationMode 'pathConstrained'.
%
%   The input must be N-by-2: [x y].  For an FIESTA Molecule use:
%       xy = Molecule(k).Results(:, 3:4);
%       S = SICstepFinder2D(xy);
%
%   Important output fields (all in localization/frame indices):
%       StepFit          N-by-2 fitted dwell coordinates
%       ChangePoints     indices i such that the step is between i and i+1
%       DwellPositions   K-by-2 mean position of each dwell
%       StepSizeStats    (K-1)-by-2 [dx dy] step vectors
%       StepMagnitudes   Euclidean lengths of the step vectors
%       StepLengths      number of points in each dwell
%       FinalSIC         SIC of the selected segmentation
%
%   joint2D detects a split from its joint x-y reduction in residual sum of
%   squares. pathConstrained uses a fitted path to suppress off-axis-only
%   excursions while retaining a full 2-D fitted mean for every dwell.

    p = inputParser;
    p.addRequired('xy', @(z) isnumeric(z) && ismatrix(z));
    p.addParameter('MinPointsInStep', 2, @(z) isnumeric(z) && isscalar(z) && z >= 2 && z == floor(z));
    p.addParameter('CorrectVariance', true, @(z) islogical(z) || (isnumeric(z) && isscalar(z)));
    p.addParameter('DecimationFactor', 1, @(z) isnumeric(z) && isscalar(z) && z >= 1 && z == floor(z));
    p.addParameter('MaxSteps', Inf, @(z) isnumeric(z) && isscalar(z) && z >= 0);
    p.addParameter('SegmentationMode', 'joint2D', @(z) ischar(z) || isstring(z));
    p.addParameter('PathModel', [], @(z) isempty(z) || isstruct(z));
    p.parse(xy, varargin{:});
    opt = p.Results;

    % A 2-by-N input is a common accidental transpose; accept it safely.
    if size(xy, 2) ~= 2 && size(xy, 1) == 2
        xy = xy.';
    end
    if size(xy, 2) ~= 2
        error('SICstepFinder2D:InputShape', 'XY must be N-by-2, with columns [x y].');
    end
    if any(~isfinite(xy(:)))
        error('SICstepFinder2D:NonFinite', 'XY contains NaN or Inf. Remove or interpolate these points first.');
    end

    mode = lower(char(opt.SegmentationMode));
    validModes = {'joint2d', 'pathconstrained'};
    if ~ismember(mode, validModes)
        error('SICstepFinder2D:SegmentationMode', ...
            'SegmentationMode must be ''joint2D'' or ''pathConstrained''.');
    end
    if strcmp(mode, 'pathconstrained') && isempty(opt.PathModel)
        error('SICstepFinder2D:PathModel', ...
            'PathModel is required when SegmentationMode is ''pathConstrained''.');
    end

    originalN = size(xy, 1);
    [xy, binEnds] = decimateMean(xy, opt.DecimationFactor);
    N = size(xy, 1);
    minPts = opt.MinPointsInStep;
    if N < minPts
        error('SICstepFinder2D:TooShort', 'XY must contain at least MinPointsInStep rows.');
    end

    % A path-constrained model uses the longitudinal coordinate for
    % changepoint selection.  This prevents lateral localization excursions
    % from creating a separate dwell, while means are still refit in x-y.
    if strcmp(mode, 'pathconstrained')
        segmentationData = projectStepFit2DPath(xy, opt.PathModel);
        likelihoodDimensions = 1;
    else
        segmentationData = xy;
        likelihoodDimensions = 2;
    end
    % Prefix sums permit an exact SSE for every trial split in O(1).
    csum = [zeros(1, size(segmentationData,2)); cumsum(segmentationData, 1)];
    csq = [0; cumsum(sum(segmentationData.^2, 2))];

    % Breaks use zero and N as sentinel boundaries.  A change at i lies
    % between rows i and i+1, exactly as in the original SICstepFinder.
    breaks = [0 N];
    if strcmp(mode, 'pathconstrained')
        % Dynamic programming finds the global minimum SSE for every number
        % of dwells, then applies SIC once across those candidate models.
        % Unlike greedy insertion, this can retain two neighboring axial
        % steps even when neither one is optimal in isolation.
        [breaks, currentSSE, currentSIC, SICHistory, SSEHistory] = ...
            globallyOptimalSIC(csum, csq, N, minPts, opt.CorrectVariance, ...
            opt.MaxSteps, likelihoodDimensions);
    else
        currentSSE = segmentationSSE(breaks, csum, csq, opt.CorrectVariance);
        currentSIC = sicGaussian(currentSSE, N, 1, likelihoodDimensions);
        SICHistory = currentSIC;
        SSEHistory = currentSSE;

        while numel(breaks) - 2 < opt.MaxSteps
            [candidate, candidateSSE] = bestSplit(breaks, csum, csq, minPts, opt.CorrectVariance);
            if isempty(candidate)
                break
            end

            proposal = sort([breaks candidate]);
            % Re-optimise nearby boundary positions after every insertion.
            [proposal, candidateSSE] = refineBoundaries(proposal, csum, csq, minPts, opt.CorrectVariance);
            nDwells = numel(proposal) - 1;
            proposalSIC = sicGaussian(candidateSSE, N, nDwells, likelihoodDimensions);

            if proposalSIC > currentSIC
                break
            end
            breaks = proposal;
            currentSSE = candidateSSE;
            currentSIC = proposalSIC;
            SICHistory(end+1, 1) = currentSIC; %#ok<AGROW>
            SSEHistory(end+1, 1) = currentSSE; %#ok<AGROW>
        end
    end

    [fitDecimated, dwellPos, stepLengthsDecimated] = makeFit(xy, breaks);
    % Expand each fitted decimated bin to its original localizations.  This
    % lets plots show raw data while retaining the intentionally coarser fit.
    fitXY = expandFit(fitDecimated, binEnds, originalN);
    changePointsDecimated = breaks(2:end-1).';
    changePoints = binEnds(changePointsDecimated);
    % Use a temporary here for compatibility with MATLAB versions that do
    % not permit colon-indexing directly on the result of another index.
    originalEnds = binEnds(breaks(2:end));
    originalBreaks = [0; originalEnds(:)];
    stepLengths = diff(originalBreaks);
    stepVectors = diff(dwellPos, 1, 1);

    % Local 2-D reduction in SSE associated with each final boundary.
    stepsChi2 = zeros(numel(changePoints), 1);
    for j = 1:numel(changePoints)
        a = breaks(j);
        b = breaks(j+1);
        c = breaks(j+2);
        stepsChi2(j) = segmentSSE(a, c, csum, csq, opt.CorrectVariance) ...
            - segmentSSE(a, b, csum, csq, opt.CorrectVariance) ...
            - segmentSSE(b, c, csum, csq, opt.CorrectVariance);
    end

    StepStatistics.StepsChi2 = stepsChi2;
    StepStatistics.NumberOfStepsFound = numel(changePoints);
    StepStatistics.StepSizeStats = stepVectors;
    StepStatistics.StepMagnitudes = hypot(stepVectors(:,1), stepVectors(:,2));
    StepStatistics.StepAngles = atan2(stepVectors(:,2), stepVectors(:,1));
    StepStatistics.StepLengths = stepLengths;
    StepStatistics.StepLengthsDecimated = stepLengthsDecimated;
    StepStatistics.StepFit = fitXY;
    StepStatistics.ChangePoints = changePoints;
    StepStatistics.ChangePointsDecimated = changePointsDecimated;
    StepStatistics.DwellPositions = dwellPos;
    StepStatistics.FinalSIC = currentSIC;
    StepStatistics.SICHistory = SICHistory;
    StepStatistics.SSEHistory = SSEHistory;
    StepStatistics.DecimatedXY = xy;
    StepStatistics.SegmentationData = segmentationData;
    StepStatistics.SegmentationMode = char(opt.SegmentationMode);
    StepStatistics.DecimationBinEnds = binEnds;
    StepStatistics.Options = opt;
end

function [bestCut, bestTotal] = bestSplit(breaks, csum, csq, minPts, corrected)
    bestCut = [];
    bestTotal = Inf;
    oldTotal = segmentationSSE(breaks, csum, csq, corrected);
    for q = 1:numel(breaks)-1
        a = breaks(q); b = breaks(q+1);
        if b - a < 2 * minPts
            continue
        end
        oldCost = segmentSSE(a, b, csum, csq, corrected);
        for cut = a + minPts : b - minPts
            total = oldTotal - oldCost ...
                + segmentSSE(a, cut, csum, csq, corrected) ...
                + segmentSSE(cut, b, csum, csq, corrected);
            if total < bestTotal
                bestTotal = total;
                bestCut = cut;
            end
        end
    end
end

function [breaks, total] = refineBoundaries(breaks, csum, csq, minPts, corrected)
    % Coordinate-descent refinement; each pass only moves a boundary if it
    % lowers the adjacent two-dwell SSE, so convergence is monotonic.
    changed = true;
    nPass = 0;
    while changed && nPass < 100
        changed = false;
        nPass = nPass + 1;
        for q = 2:numel(breaks)-1
            a = breaks(q-1); b = breaks(q); c = breaks(q+1);
            candidates = a + minPts : c - minPts;
            costs = Inf(size(candidates));
            for r = 1:numel(candidates)
                cut = candidates(r);
                costs(r) = segmentSSE(a, cut, csum, csq, corrected) + ...
                           segmentSSE(cut, c, csum, csq, corrected);
            end
            [~, ix] = min(costs);
            if candidates(ix) ~= b
                breaks(q) = candidates(ix);
                changed = true;
            end
        end
    end
    total = segmentationSSE(breaks, csum, csq, corrected);
end

function total = segmentationSSE(breaks, csum, csq, corrected)
    total = 0;
    for q = 1:numel(breaks)-1
        total = total + segmentSSE(breaks(q), breaks(q+1), csum, csq, corrected);
    end
end

function sse = segmentSSE(a, b, csum, csq, corrected)
    n = b - a;
    s = csum(b+1,:) - csum(a+1,:);
    sse = (csq(b+1) - csq(a+1)) - sum(s.^2) / n;
    % Guard against tiny negative values from floating-point cancellation.
    sse = max(sse, 0);
    if corrected && n > 1
        sse = sse * n / (n - 1);
    end
end

function value = sicGaussian(sse, N, nDwells, dimensions)
    % -2 log likelihood for an isotropic Gaussian residual plus BIC.
    % There is one fitted mean per dimension/dwell and one shared variance.
    value = dimensions * N * log(max(sse, realmin) / (dimensions * N)) + ...
        (dimensions * nDwells + 1) * log(N);
end

function [breaks, finalSSE, finalSIC, sicHistory, sseHistory] = ...
        globallyOptimalSIC(csum, csq, N, minPts, corrected, maxSteps, dimensions)
    % Find the minimum SSE segmentation for every number of dwells using
    % exact dynamic programming.  Applying SIC after this optimisation is
    % globally optimal for the supplied minimum dwell count and avoids the
    % local decisions inherent in greedy split insertion.
    maxDwells = floor(N / minPts);
    if isfinite(maxSteps)
        maxDwells = min(maxDwells, floor(maxSteps) + 1);
    end
    if maxDwells < 1
        error('SICstepFinder2D:TooShort', 'Not enough decimated data for one dwell.');
    end
    bestSSE = inf(maxDwells, N+1);
    previousBoundary = zeros(maxDwells, N+1);
    for b = minPts:N
        bestSSE(1,b+1) = segmentSSE(0, b, csum, csq, corrected);
    end
    for dwellCount = 2:maxDwells
        for b = dwellCount*minPts:N
            starts = (dwellCount-1)*minPts : b-minPts;
            candidateSSE = bestSSE(dwellCount-1, starts+1).';
            valid = isfinite(candidateSSE);
            if ~any(valid), continue, end
            starts = starts(valid);
            candidateSSE = candidateSSE(valid) + ...
                segmentSSEMany(starts, b, csum, csq, corrected);
            [bestSSE(dwellCount,b+1), location] = min(candidateSSE);
            previousBoundary(dwellCount,b+1) = starts(location);
        end
    end
    fullSSE = bestSSE(:,N+1);
    validDwells = find(isfinite(fullSSE));
    sicHistory = inf(maxDwells,1);
    for k = validDwells.'
        sicHistory(k) = sicGaussian(fullSSE(k), N, k, dimensions);
    end
    [finalSIC, dwellCount] = min(sicHistory);
    finalSSE = fullSSE(dwellCount);
    breaks = zeros(dwellCount+1,1);
    breaks(end) = N;
    b = N;
    for k = dwellCount:-1:2
        a = previousBoundary(k,b+1);
        breaks(k) = a;
        b = a;
    end
    % Keep histories compact and in increasing model-complexity order.
    sicHistory = sicHistory(validDwells);
    sseHistory = fullSSE(validDwells);
end

function sse = segmentSSEMany(a, b, csum, csq, corrected)
    a = a(:);
    n = b - a;
    sums = csum(b+1,:) - csum(a+1,:);
    sse = (csq(b+1) - csq(a+1)) - sum(sums.^2, 2) ./ n;
    sse = max(sse, 0);
    if corrected
        correction = ones(size(n));
        eligible = n > 1;
        correction(eligible) = n(eligible) ./ (n(eligible)-1);
        sse = sse .* correction;
    end
end

function [fitXY, dwellPos, lengths] = makeFit(xy, breaks)
    nDwells = numel(breaks) - 1;
    fitXY = zeros(size(xy));
    dwellPos = zeros(nDwells, 2);
    lengths = zeros(nDwells, 1);
    for q = 1:nDwells
        rows = breaks(q)+1 : breaks(q+1);
        dwellPos(q,:) = mean(xy(rows,:), 1);
        fitXY(rows,:) = repmat(dwellPos(q,:), numel(rows), 1);
        lengths(q) = numel(rows);
    end
end

function [decimatedXY, binEnds] = decimateMean(xy, factor)
    % Non-overlapping bin means reduce localization noise without discarding
    % observations.  The last bin is retained even when it is incomplete.
    N = size(xy, 1);
    binStarts = (1:factor:N).';
    binEnds = min(binStarts + factor - 1, N);
    decimatedXY = zeros(numel(binStarts), 2);
    for q = 1:numel(binStarts)
        decimatedXY(q,:) = mean(xy(binStarts(q):binEnds(q),:), 1);
    end
end

function fitXY = expandFit(fitDecimated, binEnds, originalN)
    fitXY = zeros(originalN, 2);
    startRow = 1;
    for q = 1:numel(binEnds)
        fitXY(startRow:binEnds(q),:) = repmat(fitDecimated(q,:), binEnds(q) - startRow + 1, 1);
        startRow = binEnds(q) + 1;
    end
end
