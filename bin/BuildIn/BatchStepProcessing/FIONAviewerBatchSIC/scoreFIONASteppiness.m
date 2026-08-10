function Score = scoreFIONASteppiness(data, varargin)
%SCOREFIONASTEPPINESS Assign heuristic A-F quality grades to FIONA traces.
%   Scores use the existing long/short raw and fitted traces. Regions are
%   flagged for weak step-to-noise ratio, short dwells, and (when exported)
%   EFO above a configurable multiple of the saved truncation threshold.
%   This is a review-prioritization metric, not proof that a step is real.

    p = inputParser;
    p.addParameter('DownsampleFactor', 5, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    p.addParameter('MinDwellPoints', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 1));
    p.addParameter('GoodSNR', 2, @(x) isnumeric(x) && isscalar(x) && x > 0);
    p.addParameter('EFOHighMultiplier', 2.5, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    p.addParameter('ReviewScore', 55, @(x) isnumeric(x) && isscalar(x));
    p.addParameter('GoodTracePoints', 2000, @(x) isnumeric(x) && isscalar(x) && x > 0);
    p.addParameter('GoodDirectionality', 0.90, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
    p.addParameter('MinConstantVelocityR2', 0.80, @(x) isnumeric(x) && isscalar(x) && x >= -1 && x < 1);
    p.addParameter('GoodConstantVelocityR2', 0.98, @(x) isnumeric(x) && isscalar(x) && x > -1 && x <= 1);
    p.parse(varargin{:}); opt = p.Results;
    if opt.GoodConstantVelocityR2 <= opt.MinConstantVelocityR2
        error('scoreFIONASteppiness:VelocityR2', ...
            'GoodConstantVelocityR2 must exceed MinConstantVelocityR2.');
    end
    if isempty(opt.MinDwellPoints), opt.MinDwellPoints = opt.DownsampleFactor; end
    if ~isfield(data,'trace') || size(data.trace,2) < 6
        error('scoreFIONASteppiness:Trace', 'FIONA data.trace must have six columns.');
    end
    raw = double(data.trace(:,1:2)); fit = double(data.trace(:,3:4)); N = size(raw,1);
    if isfield(data,'time') && numel(data.time) == N
        time = double(data.time(:));
    else
        time = (1:N).';
    end
    % FIONA changepoints are stored on the first sample of the new dwell.
    cp = find(any(data.trace(:,5:6) > 0,2));
    cp = cp(cp > 1);
    breaks = [0; cp(:)-1; N];
    breaks = unique(breaks, 'stable');
    boundaryMagnitude = sqrt(sum((fit(cp,:) - fit(cp-1,:)).^2,2));

    efoAvailable = isfield(data,'efo') && isfield(data,'efoThreshold') && ...
        numel(data.efo) == N && isfinite(data.efoThreshold);
    if efoAvailable
        efo = double(data.efo(:));
        if isfield(data,'efoHigh') && numel(data.efoHigh) == N
            efoHigh = logical(data.efoHigh(:));
        else
            efoHigh = highEFOMaskByMean(efo, double(data.efoThreshold), ...
                opt.EFOHighMultiplier);
        end
    else
        efoHigh = false(N,1);
    end
    Regions = repmat(emptyRegion(), numel(breaks)-1, 1);
    for j = 1:numel(Regions)
        rows = breaks(j)+1:breaks(j+1);
        residual = raw(rows,:) - fit(rows,:);
        noiseRMS = sqrt(mean(sum(residual.^2,2)));
        adjacent = [];
        if j > 1, adjacent(end+1) = boundaryMagnitude(j-1); end %#ok<AGROW>
        if j <= numel(boundaryMagnitude), adjacent(end+1) = boundaryMagnitude(j); end %#ok<AGROW>
        % An edge dwell has no step on one side. It is not evidence against
        % steppiness, so omit it from SNR grading rather than assigning zero.
        if isempty(adjacent)
            signal = NaN; snr = NaN; snrScore = 1;
        else
            signal = max(adjacent); snr = signal / max(noiseRMS, eps);
            snrScore = min(snr / opt.GoodSNR, 1);
        end
        dwellScore = min(numel(rows) / opt.MinDwellPoints, 1);
        highFraction = mean(efoHigh(rows));
        % EFO is a separate acquisition-quality warning, not a penalty that
        % turns an otherwise clear step trace into a poor-steppiness grade.
        score = 100 * (0.70*snrScore + 0.30*dwellScore);
        flags = {};
        if numel(rows) < opt.MinDwellPoints, flags{end+1} = 'short dwell'; end %#ok<AGROW>
        if ~isnan(snr) && snr < opt.GoodSNR, flags{end+1} = 'weak step/noise'; end %#ok<AGROW>
        if efoAvailable && highFraction > 0, flags{end+1} = 'high EFO'; end %#ok<AGROW>
        if isempty(flags), flagText = 'none'; else, flagText = strjoin(flags, '; '); end
        Regions(j).StartIndex = rows(1); Regions(j).EndIndex = rows(end);
        Regions(j).PointCount = numel(rows); Regions(j).NoiseRMS = noiseRMS;
        Regions(j).AdjacentStepMagnitude = signal; Regions(j).StepSNR = snr;
        Regions(j).HighEFOFraction = highFraction; Regions(j).Score = score;
        Regions(j).Grade = gradeFromScore(score); Regions(j).Flags = flagText;
        Regions(j).NeedsManualReview = score < opt.ReviewScore || highFraction > 0;
    end
    weights = [Regions.PointCount].' / max(N,1);
    validSNR = ~isnan([Regions.StepSNR].');
    if any(validSNR)
        snrWeights = weights(validSNR) / sum(weights(validSNR));
        traceSNRScore = sum(snrWeights .* min([Regions(validSNR).StepSNR].' / opt.GoodSNR, 1));
    else
        traceSNRScore = 0;
    end
    traceDwellScore = sum(weights .* min([Regions.PointCount].' / opt.MinDwellPoints, 1));
    stepEvidenceScore = 0.70 * traceSNRScore + 0.30 * traceDwellScore;
    longFit = fit(:,1);
    totalLongTravel = sum(abs(diff(longFit)));
    directionality = abs(longFit(end) - longFit(1)) / max(totalLongTravel, eps);
    directionalityScore = min(directionality / opt.GoodDirectionality, 1);
    constantVelocityR2 = linearVelocityR2(time, longFit);
    constantVelocityScore = min(max((constantVelocityR2 - opt.MinConstantVelocityR2) / ...
        (opt.GoodConstantVelocityR2 - opt.MinConstantVelocityR2), 0), 1);
    traceLengthScore = min(N / opt.GoodTracePoints, 1);
    dwellLengths = diff(breaks);
    longestDwellFraction = max(dwellLengths) / max(N,1);
    longestDwellToMedian = max(dwellLengths) / max(median(dwellLengths),1);
    % A trace must have credible steps, sufficient duration, mostly one-way
    % movement, and near-constant net progression. The latter term detects
    % extended pauses and large reversals in otherwise clean step fits.
    traceScore = 100 * (0.30 * stepEvidenceScore + 0.15 * traceLengthScore + ...
        0.15 * directionalityScore + 0.40 * constantVelocityScore);
    Trace.Score = traceScore;
    Trace.Grade = gradeFromScore(traceScore);
    Trace.NumberOfDwells = numel(Regions);
    Trace.NumberOfSteps = numel(cp);
    Trace.PointCount = N;
    Trace.StepEvidenceScore = stepEvidenceScore;
    Trace.TraceLengthScore = traceLengthScore;
    Trace.Directionality = directionality;
    Trace.DirectionalityScore = directionalityScore;
    Trace.ConstantVelocityR2 = constantVelocityR2;
    Trace.ConstantVelocityScore = constantVelocityScore;
    Trace.LongestDwellFraction = longestDwellFraction;
    Trace.LongestDwellToMedian = longestDwellToMedian;
    Trace.HighEFOFraction = mean(efoHigh);
    Trace.EFOAvailable = efoAvailable;
    Trace.NeedsManualReview = traceScore < opt.ReviewScore || Trace.HighEFOFraction > 0;
    if Trace.HighEFOFraction > 0 && traceScore < opt.ReviewScore
        Trace.Flags = 'trace progression concern; high EFO';
    elseif Trace.HighEFOFraction > 0
        Trace.Flags = 'high EFO';
    elseif traceScore < opt.ReviewScore
        Trace.Flags = 'short, non-unidirectional, or non-constant progression';
    else
        Trace.Flags = 'none';
    end
    Score.Method = 'step evidence, trace length, directionality, and constant-velocity grade; EFO reported separately';
    Score.Options = opt;
    Score.Trace = Trace;
    Score.Regions = Regions;
end

function region = emptyRegion()
    region = struct('StartIndex',[],'EndIndex',[],'PointCount',[],'NoiseRMS',[], ...
        'AdjacentStepMagnitude',[],'StepSNR',[],'HighEFOFraction',[],'Score',[], ...
        'Grade','','Flags','','NeedsManualReview',false);
end

function grade = gradeFromScore(score)
    if score >= 90, grade = 'S'; elseif score >= 80, grade = 'A'; elseif score >= 70, grade = 'B';
    elseif score >= 55, grade = 'C'; elseif score >= 40, grade = 'D'; else, grade = 'F'; end
end

function r2 = linearVelocityR2(time, position)
    valid = isfinite(time) & isfinite(position);
    if nnz(valid) < 3 || range(position(valid)) <= eps
        r2 = 0; return
    end
    coefficients = polyfit(time(valid), position(valid), 1);
    prediction = polyval(coefficients, time(valid));
    residualSum = sum((position(valid) - prediction).^2);
    totalSum = sum((position(valid) - mean(position(valid))).^2);
    r2 = max(min(1 - residualSum / max(totalSum, eps), 1), -1);
end

function high = highEFOMaskByMean(efo, threshold, multiplier)
    candidate = efo > threshold;
    padded = [false; candidate(:); false];
    starts = find(diff(padded) == 1);
    ends = find(diff(padded) == -1) - 1;
    high = false(size(candidate));
    for k = 1:numel(starts)
        rows = starts(k):ends(k);
        if mean(efo(rows)) > multiplier * threshold
            high(rows) = true;
        end
    end
end
