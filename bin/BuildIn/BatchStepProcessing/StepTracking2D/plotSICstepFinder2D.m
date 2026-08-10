function Figures = plotSICstepFinder2D(xy, StepStatistics, varargin)
%PLOTSICSTEPFINDER2D Plot a 2-D SIC step fit and its on/off-axis traces.
%
%   F = plotSICstepFinder2D(XY, S) makes two figures from an N-by-2
%   trajectory XY and the structure S returned by SICstepFinder2D:
%       1. raw xy trajectory in black with fitted dwell centres in red
%       2. on-axis and off-axis position traces with red stepwise fits
%
%   F = plotSICstepFinder2D(XY, S, 'Time', T) uses T for the horizontal
%   axis.  Otherwise frame/localization number is used.
%
%   F = plotSICstepFinder2D(..., 'OnAxis', [ux uy]) sets the on-axis
%   direction explicitly (for example, a known filament direction).  By
%   default it points from the first to the last fitted dwell centre.
%   F = plotSICstepFinder2D(..., 'PathModel', P) uses the linear or
%   quadratic centreline P returned by fitStepFit2DPath.  Its on/off values
%   are arc length and normal distance from that path.
%   F = plotSICstepFinder2D(..., 'MoleculeName', NAME) adds NAME to all
%   figure-panel titles.
%   F = plotSICstepFinder2D(..., 'MethodLabel', LABEL) changes the method
%   name in the 2-D panel title (default: SIC).
%   F = plotSICstepFinder2D(..., 'ProjectionDataStyle', STYLE) draws raw
%   projection data as a thin 'line' (default) or as 'scatter' points.
%   Horizontal grid lines on the on- and off-axis plots are spaced at 16
%   and 10 units, respectively.
%
%   Example with FIESTA data:
%       xy = Molecule(k).Results(:,3:4);
%       t  = Molecule(k).Results(:,2);
%       S = SICstepFinder2D(xy);
%       plotSICstepFinder2D(xy, S, 'Time', t);

    p = inputParser;
    p.addRequired('xy', @(z) isnumeric(z) && ismatrix(z) && size(z,2) == 2);
    p.addRequired('StepStatistics', @isstruct);
    p.addParameter('Time', [], @(z) isempty(z) || (isnumeric(z) && isvector(z)));
    p.addParameter('OnAxis', [], @(z) isempty(z) || (isnumeric(z) && numel(z) == 2));
    p.addParameter('PathModel', [], @(z) isempty(z) || isstruct(z));
    p.addParameter('MoleculeName', '', @(z) ischar(z) || isstring(z));
    p.addParameter('MethodLabel', 'SIC', @(z) ischar(z) || isstring(z));
    p.addParameter('Units', 'nm', @(z) ischar(z) || isstring(z));
    p.addParameter('Visible', 'on', @(z) ischar(z) || isstring(z));
    p.addParameter('ProjectionDataStyle', 'line', @isProjectionDataStyle);
    p.parse(xy, StepStatistics, varargin{:});
    opt = p.Results;
    moleculeName = char(opt.MoleculeName);
    projectionDataStyle = lower(char(opt.ProjectionDataStyle));
    if isempty(moleculeName)
        moleculeName = 'Molecule';
    end

    N = size(xy, 1);
    if ~isfield(StepStatistics, 'StepFit') || ~isfield(StepStatistics, 'DwellPositions')
        error('plotSICstepFinder2D:InvalidStatistics', ...
            'StepStatistics must be the output of SICstepFinder2D.');
    end
    fitXY = StepStatistics.StepFit;
    centres = StepStatistics.DwellPositions;
    if ~isequal(size(fitXY), size(xy))
        error('plotSICstepFinder2D:SizeMismatch', 'StepStatistics.StepFit must be the same size as XY.');
    end
    if isempty(opt.Time)
        time = (1:N).';
        timeLabel = 'Localization number';
    else
        time = opt.Time(:);
        if numel(time) ~= N
            error('plotSICstepFinder2D:TimeLength', 'Time must have one value per row of XY.');
        end
        timeLabel = 'Time';
    end

    if isempty(opt.PathModel) && isempty(opt.OnAxis)
        onAxis = centres(end,:) - centres(1,:);
        if norm(onAxis) < eps
            % A stationary trace has no intrinsic direction; choose PCA1.
            [~,~,V] = svd(xy - mean(xy,1), 'econ');
            onAxis = V(:,1).';
        end
    elseif isempty(opt.PathModel)
        onAxis = reshape(opt.OnAxis, 1, 2);
    end
    if isempty(opt.PathModel)
        if norm(onAxis) < eps
            error('plotSICstepFinder2D:ZeroAxis', 'OnAxis must have non-zero length.');
        end
        onAxis = onAxis / norm(onAxis);
        offAxis = [-onAxis(2), onAxis(1)];
        origin = centres(1,:);
        rawOn = (xy - origin) * onAxis.';
        rawOff = (xy - origin) * offAxis.';
        fitOn = (fitXY - origin) * onAxis.';
        fitOff = (fitXY - origin) * offAxis.';
    else
        [rawOn, rawOff] = projectStepFit2DPath(xy, opt.PathModel);
        [fitOn, fitOff] = projectStepFit2DPath(fitXY, opt.PathModel);
        onAxis = opt.PathModel.AxisU;
        offAxis = opt.PathModel.AxisV;
        origin = centres(1,:);
    end
    % Keep the first fitted dwell at zero in both projections.
    rawOn = rawOn - fitOn(1); rawOff = rawOff - fitOff(1);
    fitOn = fitOn - fitOn(1); fitOff = fitOff - fitOff(1);

    Figures.Trajectory = figure('Name', '2-D SIC step fit', 'Color', 'w', ...
        'Renderer', 'painters', 'InvertHardcopy', 'off', 'Visible', char(opt.Visible));
    ax = axes(Figures.Trajectory); hold(ax, 'on');
    plot(ax, xy(:,1), xy(:,2), 'k-', 'LineWidth', 0.75, 'DisplayName', 'Trajectory');
    plot(ax, fitXY(:,1), fitXY(:,2), '-', 'Color', [0.8 0 0], ...
        'LineWidth', 1.25, 'DisplayName', 'Step fit');
    scatter(ax, centres(:,1), centres(:,2), 33, 'r', 'filled', ...
        'MarkerEdgeColor', 'none', 'DisplayName', 'Dwell centres');
    if ~isempty(opt.PathModel)
        pathU = linspace(opt.PathModel.URange(1), opt.PathModel.URange(2), 400).';
        pathV = polyval(opt.PathModel.Coefficients, pathU);
        pathXY = opt.PathModel.Origin + pathU * opt.PathModel.AxisU + pathV * opt.PathModel.AxisV;
        plot(ax, pathXY(:,1), pathXY(:,2), '--', 'Color', [0 0.35 0.8], ...
            'LineWidth', 1, 'DisplayName', 'Path model');
    end
    axis(ax, 'equal');
    xlabel(ax, ['x (' char(opt.Units) ')']);
    ylabel(ax, ['y (' char(opt.Units) ')']);
    title(ax, sprintf('%s — 2-D %s fit: %d steps', moleculeName, char(opt.MethodLabel), StepStatistics.NumberOfStepsFound));
    applyNatureAxes(ax, false);
    styleNatureLegend(legend(ax, 'Location', 'best'));

    Figures.Projection = figure('Name', 'On/off-axis SIC step fit', 'Color', 'w', ...
        'Renderer', 'painters', 'InvertHardcopy', 'off', 'Visible', char(opt.Visible));
    % Explicit positions give the on-axis plot four times the vertical
    % space of the off-axis plot.  linkaxes below keeps time zoom/pan and
    % x-limits identical in both axes.
    left = 0.10; width = 0.86;
    axOff = axes(Figures.Projection, 'Position', [left 0.11 width 0.15]);
    axOn = axes(Figures.Projection, 'Position', [left 0.35 width 0.60]);

    hold(axOn, 'on');
    plotProjectionData(axOn, time, rawOn, projectionDataStyle, 'On-axis data');
    plot(axOn, time, fitOn, 'r-', 'LineWidth', 1.25, 'DisplayName', 'On-axis fit');
    addChangeLines(axOn, time, StepStatistics);
    ylabel(axOn, ['On-axis (' char(opt.Units) ')']);
    title(axOn, sprintf('%s — on-axis projection', moleculeName));
    applyNatureAxes(axOn, true);
    setYGridSpacing(axOn, 16);
    styleNatureLegend(legend(axOn, 'Location', 'best'));

    hold(axOff, 'on');
    plotProjectionData(axOff, time, rawOff, projectionDataStyle, 'Off-axis data');
    plot(axOff, time, fitOff, 'r-', 'LineWidth', 1.25, 'DisplayName', 'Off-axis fit');
    addChangeLines(axOff, time, StepStatistics);
    xlabel(axOff, timeLabel);
    ylabel(axOff, ['Off-axis (' char(opt.Units) ')']);
    title(axOff, sprintf('%s — off-axis projection', moleculeName));
    applyNatureAxes(axOff, true);
    setYGridSpacing(axOff, 10);
    styleNatureLegend(legend(axOff, 'Location', 'best'));
    linkaxes([axOn axOff], 'x');

    Figures.OnAxis = onAxis;
    Figures.OffAxis = offAxis;
    Figures.Origin = origin;
    Figures.PathModel = opt.PathModel;
end

function plotProjectionData(ax, time, values, style, displayName)
    switch style
        case 'line'
            plot(ax, time, values, 'k-', 'LineWidth', 0.4, 'DisplayName', displayName);
        case 'scatter'
            plot(ax, time, values, 'k.', 'MarkerSize', 5, 'DisplayName', displayName);
    end
end

function valid = isProjectionDataStyle(value)
    valid = (ischar(value) || isstring(value)) && ...
        any(strcmpi(char(value), {'line', 'scatter'}));
end

function applyNatureAxes(ax, showHorizontalGrid)
    % Compact journal-style axes suitable for vector editing and publication.
    ax.Color = 'w';
    ax.FontName = 'Arial';
    ax.FontSize = 9;
    ax.LineWidth = 0.75;
    ax.TickDir = 'out';
    ax.TickLength = [0.018 0.018];
    ax.Box = 'off';
    ax.Layer = 'top';
    ax.XColor = [0 0 0];
    ax.YColor = [0 0 0];
    ax.XGrid = 'off';
    ax.YGrid = ternary(showHorizontalGrid, 'on', 'off');
    ax.GridColor = [0.82 0.82 0.82];
    ax.GridAlpha = 0.55;
    ax.TitleFontWeight = 'normal';
    ax.TitleFontSizeMultiplier = 1;
    ax.LabelFontSizeMultiplier = 1;
end

function styleNatureLegend(handle)
    handle.FontName = 'Arial';
    handle.FontSize = 8;
    handle.Box = 'off';
    handle.Color = 'none';
end

function value = ternary(condition, whenTrue, whenFalse)
    if condition, value = whenTrue; else, value = whenFalse; end
end

function addChangeLines(ax, time, StepStatistics)
    if ~isfield(StepStatistics, 'ChangePoints') || isempty(StepStatistics.ChangePoints)
        return
    end
    cp = StepStatistics.ChangePoints(:);
    for k = 1:numel(cp)
        % Change point i is between i and i+1; place the line halfway in time.
        i = cp(k);
        if i < numel(time)
            x = (time(i) + time(i+1)) / 2;
        else
            x = time(i);
        end
        xline(ax, x, ':', 'Color', [0.7 0 0], 'HandleVisibility', 'off');
    end
end

function setYGridSpacing(ax, spacing)
    limits = ylim(ax);
    firstTick = ceil(limits(1) / spacing) * spacing;
    lastTick = floor(limits(2) / spacing) * spacing;
    ticks = firstTick:spacing:lastTick;
    if isempty(ticks)
        ticks = round(mean(limits) / spacing) * spacing;
    end
    ax.YTick = ticks;
    ax.YGrid = 'on';
end
