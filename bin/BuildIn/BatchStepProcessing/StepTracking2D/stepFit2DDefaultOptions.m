function Options = stepFit2DDefaultOptions()
%STEPFIT2DDEFAULTOPTIONS Defaults for processMoleculeStepFit2D.

    Options.Filter.Enabled = true;
    Options.Filter.MinPoints = 300;
    Options.Filter.MinNetDisplacement = 40;       % coordinate units, normally nm
    % Point-to-point path length accumulates localization noise, especially
    % in long, slow traces with extended dwells. A 0.002 cutoff retains
    % elongated, net-moving slow motors; downstream fit/score review
    % remains the safeguard against traces with poor steppiness.
    Options.Filter.MinDirectionality = 0.002;     % net displacement / path length
    Options.Filter.MinElongation = 1.2;           % PCA variance ratio

    Options.EFO.Column = 7;
    Options.EFO.ThresholdMode = 'median';
    % Threshold = median(EFO) + MADMultiplier * 1.4826 * MAD(EFO).
    % Zero means the retained trace must reach the median itself.
    Options.EFO.MADMultiplier = 0;
    Options.EFO.ConsecutivePoints = 5;
    % Reporting only: flag a remaining above-threshold EFO region only when
    % its mean EFO exceeds this multiple of the truncation threshold.
    Options.EFO.FlagMultiplier = 2.5;

    Options.Fit.MinPointsInStep = 2;
    Options.Fit.CorrectVariance = true;
    Options.Fit.DecimationFactor = 5;
    Options.Fit.MaxSteps = Inf;
    % Select 2-D dwells from motion along the fitted path, using globally
    % optimal SIC across all segmentations. This avoids off-axis blips and
    % the greedy omission of adjacent longitudinal steps.
    Options.Fit.SegmentationMode = 'pathConstrained';

    Options.Path.Degree = 2;                      % use 2 for a curved path
    Options.Alignment.Enabled = true;
    % Empty selects the SIC decimation factor in original localizations.
    % This matches FIONAviewer's AlignToMaxDx-scale local refinement.
    Options.Alignment.SearchWindow = [];
    Options.Alignment.MinDisplacement = 5;        % coordinate units
    Options.Alignment.MinPointsInDwell = 2;

    Options.Plot.InitialFits = false;
    Options.Plot.EFOAfterTruncation = true;
    Options.Plot.AlignedFits = true;
    Options.Plot.Units = 'nm';
end
