function Options = bayesStepFit2DDefaultOptions()
%BAYESSTEPFIT2DDEFAULTOPTIONS Defaults for the global Bayesian pipeline.
    Options.Filter.Enabled = true;
    Options.Filter.MinPoints = 300;
    Options.Filter.MinNetDisplacement = 40;
    Options.Filter.MinDirectionality = 0.03;
    Options.Filter.MinElongation = 1.2;
    Options.EFO.Column = 7;
    Options.EFO.ThresholdMode = 'median';
    Options.EFO.MADMultiplier = 0;
    Options.EFO.ConsecutivePoints = 5;
    Options.Bayes.MinPointsInStep = 2;
    Options.Bayes.DecimationFactor = 10;
    Options.Bayes.ExpectedDwellPoints = 30;
    Options.Bayes.PriorMeanStrength = 0.01;
    Options.Bayes.UsePathWeightedNoise = true;
    Options.Path.Degree = 1;
    Options.Alignment.Enabled = false; % preserve the Bayesian CP locations
    Options.Alignment.SearchWindow = 50;
    Options.Alignment.MinDisplacement = 5;
    Options.Alignment.MinPointsInDwell = 2;
    Options.Plot.EFOAfterTruncation = true;
    Options.Plot.FinalFits = true;
    Options.Plot.Units = 'nm';
end
