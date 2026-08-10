function Figure = plotStepFit2DEFO(time, efo, EFOInfo, moleculeName)
%PLOTSTEPFIT2DEFO Plot retained EFO values after initial high-EFO truncation.

    if nargin < 4 || isempty(moleculeName), moleculeName = 'Molecule'; end
    Figure = figure('Name', [char(moleculeName) ' EFO after truncation'], 'Color', 'w');
    ax = axes(Figure); hold(ax, 'on');
    plot(ax, time, efo, 'k.-', 'MarkerSize', 7, 'LineWidth', 0.75, 'DisplayName', 'Retained EFO');
    if isfield(EFOInfo, 'Threshold')
        thresholdLabel = 'EFO threshold';
        if isfield(EFOInfo, 'ThresholdDescription')
            thresholdLabel = EFOInfo.ThresholdDescription;
        end
        yline(ax, EFOInfo.Threshold, '--', 'Color', [0.8 0 0], 'LineWidth', 1, ...
            'DisplayName', thresholdLabel);
    end
    grid(ax, 'on'); box(ax, 'on');
    xlabel(ax, 'Time'); ylabel(ax, 'Effective frequency (EFO)');
    title(ax, sprintf('%s — EFO after truncation', char(moleculeName)));
    legend(ax, 'Location', 'best');
end
