function [Molecule, Info] = truncateHighEFOStepFit2D(Molecule, EFOOptions)
%TRUNCATEHIGHEFOSTEPFIT2D Remove initial data until EFO stays below threshold.
%   By default, threshold = median(EFO), a robust baseline estimate.  The
%   first start index with ConsecutivePoints values at or below the cutoff
%   is retained.  A nonzero MADMultiplier optionally gives a robust margin:
%   median(EFO) + MADMultiplier * 1.4826 * MAD(EFO).

    Info.Passed = false; Info.Reason = ''; Info.RemovedPoints = 0;
    if ~isfield(Molecule, 'Results') || size(Molecule.Results,2) < EFOOptions.Column
        Info.Reason = 'Missing EFO column'; return
    end
    efo = Molecule.Results(:,EFOOptions.Column);
    if any(~isfinite(efo))
        Info.Reason = 'Non-finite EFO values'; return
    end
    Info.MedianEFO = median(efo);
    Info.MADEFO = median(abs(efo - Info.MedianEFO));
    if ~isfield(EFOOptions, 'ThresholdMode') || strcmpi(EFOOptions.ThresholdMode, 'median')
        madMultiplier = 0;
        if isfield(EFOOptions, 'MADMultiplier')
            madMultiplier = EFOOptions.MADMultiplier;
        end
        Info.Threshold = Info.MedianEFO + madMultiplier * 1.4826 * Info.MADEFO;
        Info.ThresholdDescription = sprintf('median(EFO) + %.3g x robust SD', madMultiplier);
    else
        error('truncateHighEFOStepFit2D:ThresholdMode', 'Unsupported EFO threshold mode: %s', EFOOptions.ThresholdMode);
    end
    nConsecutive = EFOOptions.ConsecutivePoints;
    good = efo <= Info.Threshold;
    starts = find(conv(double(good), ones(nConsecutive,1), 'valid') == nConsecutive, 1, 'first');
    if isempty(starts)
        Info.Reason = 'EFO never remains below threshold'; return
    end
    Molecule.Results = Molecule.Results(starts:end,:);
    Info.Passed = true;
    Info.StartIndex = starts;
    Info.RemovedPoints = starts - 1;
    Info.Reason = 'Passed';
end
