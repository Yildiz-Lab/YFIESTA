function writeFIONASteppinessReports(outputFolder, Report, RegionScores)
%WRITEFIONASTEPPINESSREPORTS Write trace and region review-priority CSVs.
    valid = [Report.Succeeded];
    report = Report(valid);
    if isempty(report), return, end
    % Expand the scalar Steppiness structures into spreadsheet-friendly columns.
    traces = [report.Steppiness];
    traceTable = table({report.OutputFile}.', {traces.Grade}.', [traces.Score].', ...
        [traces.NumberOfSteps].', [traces.NumberOfDwells].', [traces.HighEFOFraction].', ...
        [traces.PointCount].', [traces.Directionality].', [traces.ConstantVelocityR2].', ...
        [traces.LongestDwellFraction].', [traces.EFOAvailable].', ...
        [traces.NeedsManualReview].', {traces.Flags}.', ...
        'VariableNames', {'OutputFile','Grade','Score','NumberOfSteps','NumberOfDwells', ...
        'HighEFOFraction','PointCount','Directionality','ConstantVelocityR2', ...
        'LongestDwellFraction','EFOAvailable','NeedsManualReview','Flags'});
    if isfield(traces,'PopulationEFO')
        screen = [traces.PopulationEFO];
        traceTable.UncappedGrade = {screen.UncappedGrade}.';
        traceTable.PopulationEFOMedian = [screen.MedianEFO].';
        traceTable.PopulationEFOReference = [screen.ReferenceEFO].';
        traceTable.PopulationEFOReferenceSource = {screen.ReferenceSource}.';
        traceTable.PopulationEFOReferenceTraceCount = [screen.ReferenceTraceCount].';
        traceTable.PopulationEFOThreshold = [screen.ThresholdEFO].';
        traceTable.PopulationEFORatio = [screen.Ratio].';
        traceTable.PopulationEFOGradeCap = {screen.GradeCap}.';
        traceTable.HighPopulationEFO = [screen.Flagged].';
        traceTable.PopulationEFOGradeCapped = [screen.GradeCapped].';
    end
    writetable(traceTable, fullfile(outputFolder, 'steppiness_trace_report.csv'));

    RegionScores = RegionScores(~cellfun(@isempty, RegionScores));
    if isempty(RegionScores), return, end
    regions = vertcat(RegionScores{:});
    if isempty(regions), return, end
    regionTable = struct2table(regions);
    writetable(regionTable, fullfile(outputFolder, 'steppiness_region_report.csv'));
end
