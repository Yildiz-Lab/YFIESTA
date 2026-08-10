function StepStatistics = BayesStepFinder2D(xy, varargin)
%BAYESSTEPFINDER2D Globally optimized Bayesian 2-D changepoint segmentation.
%   S = BAYESSTEPFINDER2D(XY) fits piecewise-constant 2-D dwell positions.
%   Unlike SICstepFinder2D, it does not greedily add the next best split.
%   It finds the maximum-a-posteriori segmentation over all allowable
%   changepoint combinations using dynamic programming.
%
%   The model has a shared 2-D localization covariance and a Gaussian prior
%   on each dwell mean.  ExpectedDwellPoints provides a geometric prior on
%   changepoints: larger values favor fewer/longer dwells.

    p = inputParser;
    p.addRequired('xy', @(z) isnumeric(z) && ismatrix(z));
    p.addParameter('MinPointsInStep', 2, @(z) isnumeric(z) && isscalar(z) && z >= 2 && z == floor(z));
    p.addParameter('DecimationFactor', 10, @(z) isnumeric(z) && isscalar(z) && z >= 1 && z == floor(z));
    p.addParameter('ExpectedDwellPoints', 30, @(z) isnumeric(z) && isscalar(z) && z > 1);
    p.addParameter('PriorMeanStrength', 0.01, @(z) isnumeric(z) && isscalar(z) && z > 0);
    p.addParameter('NoiseCovariance', [], @(z) isempty(z) || (isnumeric(z) && isequal(size(z),[2 2])));
    p.parse(xy, varargin{:}); opt = p.Results;

    if size(xy,2) ~= 2 && size(xy,1) == 2, xy = xy.'; end
    if size(xy,2) ~= 2 || any(~isfinite(xy(:)))
        error('BayesStepFinder2D:Input', 'XY must be a finite N-by-2 array.');
    end
    originalN = size(xy,1);
    [data, binEnds] = decimateMean(xy, opt.DecimationFactor);
    N = size(data,1); minPts = opt.MinPointsInStep;
    if N < minPts
        error('BayesStepFinder2D:TooShort', 'Too few decimated points for MinPointsInStep.');
    end

    if isempty(opt.NoiseCovariance)
        covariance = robustCovariance(data);
    else
        covariance = (opt.NoiseCovariance + opt.NoiseCovariance.') / 2;
    end
    [cholCov,pd] = chol(covariance);
    if pd ~= 0
        error('BayesStepFinder2D:Covariance', 'NoiseCovariance must be positive definite.');
    end
    inverseCovariance = cholCov \ (cholCov' \ eye(2));
    logDetCovariance = 2 * sum(log(diag(cholCov)));
    priorMean = mean(data,1);
    kappa0 = opt.PriorMeanStrength;
    hazard = 1 / opt.ExpectedDwellPoints;
    changePenalty = log((1-hazard) / hazard);

    sums = [zeros(1,2); cumsum(data,1)];
    outer = zeros(N+1,3); % cumulative [x^2 xy y^2]
    outer(2:end,:) = cumsum([data(:,1).^2 data(:,1).*data(:,2) data(:,2).^2],1);
    score = -Inf(N,1); previous = zeros(N,1);

    % Dynamic programming: every segmentation is considered, so this is a
    % global MAP solution rather than sequential greedy split insertion.
    for endIndex = minPts:N
        bestScore = -Inf; bestStart = 0;
        for startIndex = 1:endIndex-minPts+1
            segmentScore = logSegmentMarginal(startIndex,endIndex,sums,outer, ...
                inverseCovariance,logDetCovariance,priorMean,kappa0);
            if startIndex == 1
                candidate = segmentScore;
            elseif isfinite(score(startIndex-1))
                candidate = score(startIndex-1) + segmentScore - changePenalty;
            else
                continue
            end
            if candidate > bestScore
                bestScore = candidate; bestStart = startIndex;
            end
        end
        score(endIndex) = bestScore; previous(endIndex) = bestStart;
    end
    if ~isfinite(score(N))
        error('BayesStepFinder2D:Segmentation', 'No valid segmentation was found.');
    end

    breaks = N;
    endIndex = N;
    while previous(endIndex) > 1
        startIndex = previous(endIndex);
        breaks = [startIndex-1; breaks]; %#ok<AGROW>
        endIndex = startIndex-1;
    end
    breaks = [0; breaks];
    [fitDecimated, centres, decimatedLengths] = makeFit(data, breaks);
    fitXY = expandFit(fitDecimated, binEnds, originalN);
    changePointsDecimated = breaks(2:end-1);
    changePoints = binEnds(changePointsDecimated);
    originalEnds = binEnds(breaks(2:end));
    originalBreaks = [0; originalEnds(:)];
    vectors = diff(centres,1,1);

    StepStatistics.Method = 'Bayesian global MAP changepoint segmentation';
    StepStatistics.NumberOfStepsFound = numel(changePoints);
    StepStatistics.ChangePoints = changePoints;
    StepStatistics.ChangePointsDecimated = changePointsDecimated;
    StepStatistics.StepFit = fitXY;
    StepStatistics.DwellPositions = centres;
    StepStatistics.StepSizeStats = vectors;
    StepStatistics.StepMagnitudes = hypot(vectors(:,1),vectors(:,2));
    StepStatistics.StepAngles = atan2(vectors(:,2),vectors(:,1));
    StepStatistics.StepLengths = diff(originalBreaks);
    StepStatistics.StepLengthsDecimated = decimatedLengths;
    StepStatistics.StepsChi2 = localBayesEvidence(breaks,sums,outer,inverseCovariance, ...
        logDetCovariance,priorMean,kappa0);
    StepStatistics.FinalLogPosterior = score(N);
    StepStatistics.NoiseCovariance = covariance;
    StepStatistics.DecimatedXY = data;
    StepStatistics.DecimationBinEnds = binEnds;
    StepStatistics.Options = opt;
end

function value = logSegmentMarginal(a,b,sums,outer,inverseCov,logDetCov,mu0,kappa0)
    n = b-a+1;
    sumX = sums(b+1,:) - sums(a,:);
    meanX = sumX/n;
    q = outer(b+1,:) - outer(a,:);
    scatter = [q(1) q(2); q(2) q(3)] - n*(meanX.'*meanX);
    delta = meanX - mu0;
    kappaN = kappa0 + n;
    quadratic = trace(inverseCov*scatter) + (kappa0*n/kappaN) * (delta*inverseCov*delta.');
    value = -n*log(2*pi) - 0.5*n*logDetCov + log(kappa0/kappaN) - 0.5*quadratic;
end

function covariance = robustCovariance(data)
    differences = diff(data,1,1);
    if isempty(differences), differences = data - mean(data,1); end
    scale = 1.4826 * median(abs(differences - median(differences,1)),1) / sqrt(2);
    scale = max(scale, eps(max(abs(data(:)))));
    covariance = diag(scale.^2);
end

function evidence = localBayesEvidence(breaks,sums,outer,inverseCov,logDetCov,mu0,kappa0)
    evidence = zeros(numel(breaks)-2,1);
    for j = 1:numel(evidence)
        a = breaks(j)+1; b = breaks(j+1); c = breaks(j+2);
        evidence(j) = logSegmentMarginal(a,b,sums,outer,inverseCov,logDetCov,mu0,kappa0) + ...
            logSegmentMarginal(b+1,c,sums,outer,inverseCov,logDetCov,mu0,kappa0) - ...
            logSegmentMarginal(a,c,sums,outer,inverseCov,logDetCov,mu0,kappa0);
    end
end

function [decimatedXY, binEnds] = decimateMean(xy, factor)
    N = size(xy,1); starts = (1:factor:N).'; binEnds = min(starts+factor-1,N);
    decimatedXY = zeros(numel(starts),2);
    for j = 1:numel(starts), decimatedXY(j,:) = mean(xy(starts(j):binEnds(j),:),1); end
end

function [fitXY,centres,lengths] = makeFit(xy,breaks)
    nDwells = numel(breaks)-1; fitXY = zeros(size(xy)); centres = zeros(nDwells,2); lengths = zeros(nDwells,1);
    for j = 1:nDwells
        rows = breaks(j)+1:breaks(j+1); centres(j,:) = mean(xy(rows,:),1);
        fitXY(rows,:) = repmat(centres(j,:),numel(rows),1); lengths(j) = numel(rows);
    end
end

function fitXY = expandFit(fitDecimated,binEnds,originalN)
    fitXY = zeros(originalN,2); startRow = 1;
    for j = 1:numel(binEnds)
        fitXY(startRow:binEnds(j),:) = repmat(fitDecimated(j,:),binEnds(j)-startRow+1,1);
        startRow = binEnds(j)+1;
    end
end
