function Figure = plotFIONASteppinessFlags(source, varargin)
%PLOTFIONASTEPPINESSFLAGS Plot FIONA traces with manual-review regions.
%   Black = no flag, red = weak/short stepping region, blue = high-EFO
%   region. Blue has precedence when a region has both types of flag.
%
%   plotFIONASteppinessFlags(DATA) accepts a loaded FIONA data structure.
%   plotFIONASteppinessFlags(FILE) loads a *_fiona.mat file directly.
%   plotFIONASteppinessFlags() opens a file chooser. Both the 2-D SIC
%   pre-screen field (data.Steppiness2DPreScreen) and the completed
%   FIONAviewer field (data.Steppiness) are supported.

    if nargin < 1 || isempty(source)
        [filename, folder] = uigetfile({'*_fiona.mat', 'FIONA trace MAT files'}, ...
            'Select a FIONA trace for steppiness review');
        if isequal(filename, 0)
            Figure = [];
            return
        end
        source = fullfile(folder, filename);
    end
    [data, score, defaultName] = resolveFIONASource(source);

    p = inputParser;
    p.addParameter('Name', 'FIONA steppiness review', @(x) ischar(x) || isstring(x));
    p.addParameter('Visible', 'on', @(x) ischar(x) || isstring(x));
    p.parse(varargin{:}); opt = p.Results;
    if ismember('Name', p.UsingDefaults) && ~isempty(defaultName)
        opt.Name = defaultName;
    end
    if ~isfield(data,'trace') || size(data.trace,2) < 6
        error('plotFIONASteppinessFlags:Data', ...
            'The selected FIONA data must contain data.trace with six columns.');
    end
    N = size(data.trace,1);
    if isfield(data,'time') && numel(data.time) == N
        time = double(data.time(:)); xLabel = 'time';
    else
        time = (1:N).'; xLabel = 'frame';
    end
    % 0 = normal (black), 1 = steppiness review (red), 2 = high EFO (blue).
    class = zeros(N,1);
    regions = score.Regions;
    reviewScore = score.Options.ReviewScore;
    for k = 1:numel(regions)
        if regions(k).Score < reviewScore
            class(regions(k).StartIndex:regions(k).EndIndex) = 1;
        end
    end
    if isfield(data,'efoHigh') && numel(data.efoHigh) == N
        class(logical(data.efoHigh(:))) = 2;
    else
        for k = 1:numel(regions)
            if regions(k).HighEFOFraction > 0
                class(regions(k).StartIndex:regions(k).EndIndex) = 2;
            end
        end
    end

    Figure = figure('Name', char(opt.Name), 'Color', 'w', 'Visible', char(opt.Visible), ...
        'Renderer', 'painters', 'Position', [100 100 900 600]);
    % Match the projection layout used for the 2-D SIC figures: long axis
    % receives four-fifths of the height and short axis one-fifth.
    layout = tiledlayout(Figure, 5, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    axesHandles = gobjects(2,1);
    labels = {'Long axis', 'Short axis'};
    for axisIndex = 1:2
        if axisIndex == 1
            axesHandles(axisIndex) = nexttile(layout, [4 1]);
        else
            axesHandles(axisIndex) = nexttile(layout);
        end
        raw = double(data.trace(:,axisIndex));
        fit = double(data.trace(:,axisIndex+2));
        hold on
        plot(time, raw, '-', 'Color', [0 0 0], 'LineWidth', 0.7, 'HandleVisibility', 'off');
        plot(time, fit, '-', 'Color', [0 0 0], 'LineWidth', 1.5, 'HandleVisibility', 'off');
        plotFlaggedSegments(time, raw, fit, class, 1, [0.85 0.10 0.10]);
        plotFlaggedSegments(time, raw, fit, class, 2, [0.10 0.35 0.90]);
        ylabel(labels{axisIndex}); grid on; box off
        if axisIndex == 1, title(char(opt.Name), 'Interpreter', 'none'); end
        if axisIndex == 2, xlabel(xLabel); end
    end
    linkaxes(axesHandles, 'x');
    % Keep the legend intentionally limited to the two review flags.
    blueKey = plot(axesHandles(1), nan, nan, '-', 'Color', [0.10 0.35 0.90], ...
        'LineWidth', 2.2, 'DisplayName', 'high EFO (mean > 2.5 x threshold)');
    redKey = plot(axesHandles(1), nan, nan, '-', 'Color', [0.85 0.10 0.10], ...
        'LineWidth', 2.2, 'DisplayName', 'short dwell or weak step');
    legend(axesHandles(1), [blueKey redKey], 'Location', 'best');
end

function [data, score, displayName] = resolveFIONASource(source)
    displayName = '';
    if ischar(source) || isstring(source)
        filename = char(source);
        if ~isfile(filename)
            error('plotFIONASteppinessFlags:File', 'FIONA file does not exist: %s', filename);
        end
        loaded = load(filename);
        if isfield(loaded, 'data') && isstruct(loaded.data)
            data = loaded.data;
        elseif isfield(loaded, 'Data') && isstruct(loaded.Data)
            data = loaded.Data;
        else
            names = fieldnames(loaded);
            data = [];
            for k = 1:numel(names)
                value = loaded.(names{k});
                if isstruct(value) && isfield(value, 'trace')
                    data = value;
                    break
                end
            end
            if isempty(data)
                error('plotFIONASteppinessFlags:FileFormat', ...
                    'No FIONA data structure was found in %s.', filename);
            end
        end
        [~, stem] = fileparts(filename);
        displayName = [stem ' - steppiness review'];
    elseif isstruct(source)
        data = source;
    else
        error('plotFIONASteppinessFlags:Input', ...
            'Input must be a FIONA data structure or a *_fiona.mat filename.');
    end

    if isfield(data, 'Steppiness') && isstruct(data.Steppiness)
        score = data.Steppiness;
    elseif isfield(data, 'Steppiness2DPreScreen') && isstruct(data.Steppiness2DPreScreen)
        score = data.Steppiness2DPreScreen;
    else
        % Permit older or manually edited FIONA traces to be reviewed too.
        % This does not alter the MAT file; it only computes a temporary
        % score for the current figure.
        score = scoreFIONASteppiness(data);
    end
    if ~isfield(score, 'Regions') || ~isfield(score, 'Options') || ...
            ~isfield(score.Options, 'ReviewScore')
        error('plotFIONASteppinessFlags:Score', ...
            'The selected FIONA trace does not contain a usable steppiness score.');
    end
end

function plotFlaggedSegments(time, raw, fit, class, code, color)
    selected = class == code;
    padded = [false; selected; false];
    starts = find(diff(padded) == 1);
    ends = find(diff(padded) == -1) - 1;
    for k = 1:numel(starts)
        rows = starts(k):ends(k);
        plot(time(rows), raw(rows), '-', 'Color', color, 'LineWidth', 0.9, ...
            'HandleVisibility', 'off');
        plot(time(rows), fit(rows), '-', 'Color', color, 'LineWidth', 2.2, ...
            'HandleVisibility', 'off');
    end
end
