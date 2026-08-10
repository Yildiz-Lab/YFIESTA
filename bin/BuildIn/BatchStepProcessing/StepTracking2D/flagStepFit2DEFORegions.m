function EFOFlags = flagStepFit2DEFORegions(time, efo, threshold, flagMultiplier)
%FLAGSTEPFIT2DEFOREGIONS Flag unusually high contiguous EFO regions.
%   The threshold is the robust median-based EFO cutoff used during the
%   StepFit2D preparation stage. A contiguous candidate region must be above
%   that threshold and have mean EFO > FLAGMULTIPLIER*threshold (default 2.5)
%   to become a review flag. These flags never alter a completed step fit.

    if nargin < 4 || isempty(flagMultiplier), flagMultiplier = 2.5; end
    time = double(time(:)); efo = double(efo(:));
    if numel(time) ~= numel(efo)
        error('flagStepFit2DEFORegions:Length', 'time and efo must have equal length.');
    end
    candidateHigh = efo > threshold;
    padded = [false; candidateHigh; false];
    starts = find(diff(padded) == 1);
    ends = find(diff(padded) == -1) - 1;
    Regions = repmat(struct('StartIndex',[],'EndIndex',[],'StartTime',[], ...
        'EndTime',[],'PointCount',[],'MeanEFO',[],'MaxEFO',[]), numel(starts), 1);
    high = false(size(efo));
    kept = 0;
    for k = 1:numel(starts)
        rows = starts(k):ends(k);
        meanEFO = mean(efo(rows));
        if meanEFO <= flagMultiplier * threshold, continue, end
        kept = kept + 1;
        high(rows) = true;
        Regions(kept).StartIndex = rows(1); Regions(kept).EndIndex = rows(end);
        Regions(kept).StartTime = time(rows(1)); Regions(kept).EndTime = time(rows(end));
        Regions(kept).PointCount = numel(rows);
        Regions(kept).MeanEFO = meanEFO; Regions(kept).MaxEFO = max(efo(rows));
    end
    Regions = Regions(1:kept);
    EFOFlags.Threshold = threshold;
    EFOFlags.FlagMultiplier = flagMultiplier;
    EFOFlags.FlagMeanThreshold = flagMultiplier * threshold;
    EFOFlags.HighMask = high;
    EFOFlags.HighFraction = mean(high);
    EFOFlags.Regions = Regions;
end
