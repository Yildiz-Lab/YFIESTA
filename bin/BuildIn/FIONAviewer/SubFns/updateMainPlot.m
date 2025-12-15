function app = updateMainPlot(app)
    % Check that data exists
    
    hold(app.AxMain, 'on')
    cla(app.AxMain);
    if ~isempty(app.Data.PSD1Data_Long) && isfield(app.Data,'t')
        % plot data
        if app.Data.FilteredFlag %show filtered if toggled (won't save filtered though, pulls from raw data for that one)
            if strcmp(app.DisplayPanelFields.DropLineStyle.Value, 'Lines')
                hLong = plot(app.AxMain, app.Data.t_Filt, app.Data.PSD1Data_Long_Filt, 'Tag', 'hLongAxis');
                hShort = plot(app.AxMain, app.Data.t_Filt, app.Data.PSD1Data_Short_Filt, 'Tag', 'hShortAxis');
            elseif strcmp(app.DisplayPanelFields.DropLineStyle.Value, 'Dots')
                hLong = plot(app.AxMain, app.Data.t_Filt, app.Data.PSD1Data_Long_Filt,'.', 'Tag', 'hLongAxis');
                hShort = plot(app.AxMain, app.Data.t_Filt, app.Data.PSD1Data_Short_Filt,'.', 'Tag', 'hShortAxis');
            end
        else
            if strcmp(app.DisplayPanelFields.DropLineStyle.Value, 'Lines')
                hLong = plot(app.AxMain, app.Data.t, app.Data.PSD1Data_Long, 'Tag', 'hLongAxis');
                hShort = plot(app.AxMain, app.Data.t, app.Data.PSD1Data_Short, 'Tag', 'hShortAxis');
            elseif strcmp(app.DisplayPanelFields.DropLineStyle.Value, 'Dots')
                hLong = plot(app.AxMain, app.Data.t, app.Data.PSD1Data_Long,'.', 'Tag', 'hLongAxis');
                hShort = plot(app.AxMain, app.Data.t, app.Data.PSD1Data_Short,'.', 'Tag', 'hShortAxis');
            end
        end
        
        % make sure all lines are not pickable for point-click
        % functionality
        hLong.PickableParts = 'none';
        hShort.PickableParts = 'none';
        
        % Now plot steps
        if ~isempty(app.Data.stepVector)
            hFit = plot(app.AxMain, app.Data.t, app.Data.stepVector, 'Color', [80,200,120]/255, 'LineWidth', 1.5, 'Tag', 'hFitLine');
            hFit.PickableParts = 'none';
        end
        if ~isempty(app.Data.shortstepVector)
            hFitShort = plot(app.AxMain, app.Data.t, app.Data.shortstepVector, 'Color', [50,50,50]/255, 'LineWidth', 1.5, 'Tag', 'hFitShortLine');
            hFitShort.PickableParts = 'none';
        end
        
        
        if app.time_bool
            xlabel(app.AxMain,'time (s)');
        else
            xlabel(app.AxMain,'frame');
        end
        ylabel(app.AxMain,'position (nm)');
        title(app.AxMain,'PSD1 Data');
        
        % Show grid option or not
        if app.DisplayPanelFields.ChkGrid.Value
            % grid(app.AxMain,'on');
            dx = app.DisplayPanelFields.GridSpacing.Value;
            if ~isempty(app.Data.PSD1Data_Long)
                mn = min(app.Data.PSD1Data_Short,app.Data.PSD1Data_Long,'omitnan');
                mx = max(app.Data.PSD1Data_Short,app.Data.PSD1Data_Long,'omitnan');
                mn = min([app.AxMain.YLim(1),mn])-dx; mx = max([app.AxMain.YLim(2),mx])+dx;
            else
                mn = app.AxMain.YLim(1); mx = app.AxMain.YLim(2);
            end
            set(app.AxMain,'YTick', sort([0:-dx:mn, dx:dx:mx]))
            app.AxMain.YGrid = 'on'; app.AxMain.XGrid = 'off';
        else
            grid(app.AxMain,'off');
        end

    else
        % Clear axes if no data
        cla(app.AxMain);
        title(app.AxMain,'No data loaded');
        
    end
    hold(app.AxMain, 'off')

    title(app.AxMain, app.FileName.Text, 'Interpreter', 'None')

end