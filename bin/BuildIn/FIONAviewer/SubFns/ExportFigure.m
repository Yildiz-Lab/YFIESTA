function ExportFigure(app)

    % Set everything to pixels
    app.AxMain.Parent.Units = 'pixels';
    app.AxMain.Units = 'pixels';

    % Create new figure
    figWidth = 600; % or any desired figure width
    figHeight = 400; % desired figure height
    figExport = figure('Name','Exported trace','Color','w', ...
                       'Units','pixels', 'Position',[100 100 figWidth figHeight]);

    % Create axes that fills the figure
    axNew = axes(figExport, 'Units','normalized', 'Position',[0.1 0.1 0.85 0.85]); 
    % normalized position with small margins for labels

    hold(axNew, 'on');

    % Replot all lines
    lines = findobj(app.AxMain, 'Type', 'Line');
    for k = 1:numel(lines)
        plot(axNew, lines(k).XData, lines(k).YData, ...
            'Color', lines(k).Color, ...
            'LineStyle', lines(k).LineStyle, ...
            'Marker', lines(k).Marker, ...
            'LineWidth', lines(k).LineWidth, ...
            'MarkerSize', lines(k).MarkerSize);
    end

    % Copy axes properties
    axNew.XLim = app.AxMain.XLim;
    axNew.YLim = app.AxMain.YLim;
    axNew.XScale = app.AxMain.XScale;
    axNew.YScale = app.AxMain.YScale;
    axNew.Box = app.AxMain.Box;
    axNew.Color = app.AxMain.Color;
    axNew.XTick = app.AxMain.XTick;
    axNew.YTick = app.AxMain.YTick;
    axNew.XGrid = app.AxMain.XGrid;
    axNew.YGrid = app.AxMain.YGrid;
    axNew.XMinorGrid = app.AxMain.XMinorGrid;
    axNew.YMinorGrid = app.AxMain.YMinorGrid;
    axNew.TickDir = app.AxMain.TickDir;
    axNew.FontSize = app.AxMain.FontSize;
    axNew.FontName = app.AxMain.FontName;

    % Copy labels and title
    axNew.XLabel.String = app.AxMain.XLabel.String;
    axNew.YLabel.String = app.AxMain.YLabel.String;
    axNew.Title.String  = app.AxMain.Title.String;
    axNew.XLabel.FontSize = app.AxMain.XLabel.FontSize;
    axNew.YLabel.FontSize = app.AxMain.YLabel.FontSize;
    axNew.Title.FontSize  = app.AxMain.Title.FontSize;
    axNew.Title.Interpreter = 'none';

    hold(axNew, 'off');

    % Automatically save to figure
    [pathname,filename,~] = fileparts(app.LeftPanelFields.StepsFilename.Text);
    saveas(figExport , fullfile(pathname,strcat(filename,'.fig')));
    close(figExport);

    disp(strcat("Figure exported to ", fullfile(pathname,strcat(filename,'.fig'))))

    figure(app.fig);