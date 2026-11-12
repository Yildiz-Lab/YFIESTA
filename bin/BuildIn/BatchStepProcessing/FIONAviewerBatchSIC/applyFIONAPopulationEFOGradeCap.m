function Trace = applyFIONAPopulationEFOGradeCap(Trace, Screen, index)
%APPLYFIONAPOPULATIONEFOGRADECAP Cap the grade without changing fit or score.
    grades = {'S','A','B','C','D','F'};
    if isfield(Trace,'PopulationEFO')
        previous = Trace.PopulationEFO;
        Trace.Grade = previous.UncappedGrade;
        Trace.Flags = previous.UncappedFlags;
        Trace.NeedsManualReview = previous.UncappedNeedsManualReview;
    end
    info = struct('UncappedGrade',Trace.Grade, 'UncappedFlags',Trace.Flags, ...
        'UncappedNeedsManualReview',Trace.NeedsManualReview, ...
        'MedianEFO',Screen.MedianEFO(index), 'ReferenceEFO',Screen.ReferenceEFO, ...
        'ReferenceTraceCount',Screen.ReferenceTraceCount, 'ReferenceSource',Screen.ReferenceSource, ...
        'Multiplier',Screen.Multiplier, 'ThresholdEFO',Screen.ThresholdEFO, ...
        'Ratio',Screen.Ratio(index), 'GradeCap',Screen.GradeCap, ...
        'Flagged',Screen.Flagged(index), 'GradeCapped',false);
    if info.Flagged
        original = find(strcmp(grades, Trace.Grade),1);
        cap = find(strcmp(grades, info.GradeCap),1);
        if ~isempty(original) && original < cap
            Trace.Grade = grades{cap};
            info.GradeCapped = true;
        end
        Trace.NeedsManualReview = true;
        flag = sprintf('high population EFO (%.2fx reference)',info.Ratio);
        if isempty(Trace.Flags) || strcmp(Trace.Flags,'none')
            Trace.Flags = flag;
        else
            Trace.Flags = [Trace.Flags '; ' flag];
        end
    end
    Trace.PopulationEFO = info;
end
