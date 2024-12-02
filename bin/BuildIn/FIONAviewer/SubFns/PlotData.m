function handles = PlotData(hObject, handles, keepLimits, neighbors)
% Plots the trap data with all the current settings.
% Generates the following variables in handles: currentPlotT, 
% currentPlotPSD_Long, currentPlotPSD_Short, currentPlotTrap_Long,
% and currentPlotTrap_Short
% If keepLimits = 1, attempts to keep previously selected Y-limits
% If keepLimits = 2, keeps the limits exactly the same
% Otherwise, resets the Y-limits to default values
%
% Modified from Tom Bilyard's files PlotMainData.m and SaveYScale.m 
% by Vladislav Belyy
% Last updated on 11/18/2011

%JS Edit 2022/09/27 Add neighbors into plot


%method=1 for PS, =0 for Data display

%% Generate data to be plotted

doNotCenter = 0; 
allowTrapPosDisplay = 0; % only allowed if position signals are in nm and 
 % rotated; otherwise, displaying trap position makes no sense

if keepLimits % keep the same limits despite changing units
  
    oldXlim=get(handles.axes1,'Xlim'); % save old X limit
    oldYlim=get(handles.axes1,'Ylim'); % save old Y limit
    
    minPos = min(handles.currentPlotPSD_Long);
    maxPos = max(handles.currentPlotPSD_Long);
    
    minT = min(handles.currentPlotT);
    maxT = max(handles.currentPlotT);
end

    

% set to always plot raw data0
% plot filtered data if you find any
%can do something similar to allow centering support for steps

%can also not use _2 and just use filtered data right off

%should go make sure times are copasetic ion filter m-file

% if handles.display.filtered % Filtered or decimated trace
% 
%     handles.currentPlotPSD_Long_2 = handles.PSD1Data_Long_Filt;
%     handles.currentPlotPSD_Short_2 = handles.PSD1Data_Short_Filt;
%     
% 
%     %is this handle made in the filter setup 
%     handles.currentPlotT_2 = handles.t;
% end
    
handles.currentPlotPSD_Long = handles.PSD1Data_Long;
handles.currentPlotPSD_Short = handles.PSD1Data_Short;

handles.currentPlotTrap_Long = handles.trapPosLong;
handles.currentPlotTrap_Short = handles.trapPosShort;

handles.currentPlotT = handles.t;



handles.currentYlabel='position (nm)';
        



%% Center traces, if necessary
%might not work with steps trace, but put filtered bits in there for now
if ~doNotCenter 
    if (handles.display.center == 1) % median centering
        medianValueX = median(RemoveNaNs(handles.currentPlotPSD_Long));
        medianValueY = median(RemoveNaNs(handles.currentPlotPSD_Short));

        handles.currentPlotPSD_Long = ...
            handles.currentPlotPSD_Long - medianValueX;

        handles.currentPlotPSD_Short = ...
            handles.currentPlotPSD_Short - medianValueY;
%         if handles.display.filtered        
%             handles.currentPlotPSD_Long_2 = ...
%                 handles.currentPlotPSD_Long_2 - medianValueX;
% 
%             handles.currentPlotPSD_Short_2 = ...
%                 handles.currentPlotPSD_Short_2 - medianValueX;
%         end
%             
    
    elseif (handles.display.center == 2) % user value centering
        centerValueX = str2double(get(handles.XOffset,'String'));
        centerValueY = str2double(get(handles.YOffset,'String'));

        handles.currentPlotPSD_Long = ...
            handles.currentPlotPSD_Long - centerValueX;

        handles.currentPlotPSD_Short = ...
            handles.currentPlotPSD_Short - centerValueY;
%         if handles.display.filtered
%             handles.currentPlotPSD_Long_2 = ...
%                 handles.currentPlotPSD_Long_2 - centerValueX;
% 
%             handles.currentPlotPSD_Short_2 = ...
%                 handles.currentPlotPSD_Short_2 - centerValueY;
%         end
    end
end



%% Adjust display limits

if keepLimits == 1 % attempt to keep the same limits despite changing units
    
    minPosNew = min(real(handles.currentPlotPSD_Long));
    maxPosNew = max(real(handles.currentPlotPSD_Long));
    
    minTNew = min(handles.currentPlotT);
    maxTNew = max(handles.currentPlotT);
    
    % linear transformation; newLim = m*(oldLim) + b
    % first, find coeffs mx and bx for position and mt and bt for time
    mx =  (minPosNew - maxPosNew) / (minPos - maxPos);
    bx = minPosNew - mx*minPos;
    
    mt =  (minTNew - maxTNew) / (minT - maxT);
    bt = minTNew - mx*minT;
    
    % Transform old limits into new limits:
    newXlim = sort(real(mt*oldXlim + bt));
    newYlim = sort(real(mx*oldYlim + bx));

elseif keepLimits == 2 % Keep limits strictly the same
    
    newXlim = oldXlim;
    newYlim = oldYlim;

end
    
    




%% Actual plotting

% Read plotting parameters
filename = get(handles.FileName,'string');
dispLong = get(handles.ShowX,'value'); % Plot long axis?
dispShort = get(handles.ShowY,'value'); % plot short axis?
if handles.display.filtered
    dispLongFilt = get(handles.ShowXFilt,'value'); %plot long filtered axis
    dispShortFilt = get(handles.ShowYFilt,'value'); %plot short filtered axis
else
    dispLongFilt = 0;
    dispShortFilt = 0;
end

plotType=get(handles.LinesOrDots,'string');
val=get(handles.LinesOrDots,'value');

switch plotType{val};
case 'Dots' % User selects dots
   handles.currentstyleX='b.';
   handles.currentstyleY='r.';
   handles.currentstyleXFilt = 'g.';
   handles.currentstyleYFilt = 'g.';
case 'Lines' % User selects lines
   handles.currentstyleX = 'b-';
   handles.currentstyleY = 'r-';
   handles.currentstyleXFilt = 'g-';
   handles.currentstyleYFilt = 'g-';
end

axes(handles.axes1) %#ok<MAXES>

plot(0,0)

hold on

if handles.display.gridLines % display gridlines
    
    % Determine user-specified gridline spacing
    DivString=get(handles.GridDiv,'string');
    Loc=strfind(DivString,' ');

    if size(Loc,1)==0 % No units specified
        Div = str2double(DivString);
    else
        Div = str2double(DivString(1:strfind(DivString,' ')-1));
    end

    %Ylims=get(handles.axes1,'Ylim'); % get Y-limits
    
    
    % Generate ticks
    minY = min(min(real(handles.currentPlotPSD_Long)), ...
        min(real(handles.currentPlotPSD_Short)));
    maxY = max(max(real(handles.currentPlotPSD_Long)), ...
        max(real(handles.currentPlotPSD_Short)));
    
    % If trap position is displayed, include minimum and maximum of that
%     if handles.display.trapPos && allowTrapPosDisplay
%         minY = min(minY, min(handles.currentPlotTrap_Long));
%         maxY = max(maxY, max(handles.currentPlotTrap_Long));
%     end
    
    Ticks=floor(minY/Div)*Div:Div:ceil(maxY/Div)*Div;
    
    set(handles.axes1,'YGrid','on')
    set(handles.axes1,'YTick',Ticks)

else % Do not display gridlines
    set(handles.axes1,'YGrid','off')

end

FilteredPlotOffset = str2double(get(handles.FiltOffsetDistance,'String'));


if dispShort == 1 % Display short axis

    plot(real(handles.currentPlotPSD_Short), ...
        handles.currentstyleY)

end

if dispLong == 1 % Display long axis

    plot(real(handles.currentPlotPSD_Long), ...
        handles.currentstyleX)

end

if dispLongFilt == 1 && handles.FilteredFlag == 1   %Display filtered long axis
    
    plot(real(handles.PSD1Data_Long_Filt+FilteredPlotOffset), ...
        handles.currentstyleXFilt)
    
end

if dispShortFilt == 1 && handles.FilteredFlag == 1  %Display filtered short axis
    
    plot(real(handles.PSD1Data_Short_Filt+FilteredPlotOffset), ...
        handles.currentstyleYFilt)
    
end

xlabel(handles.currentXlabel);
ylabel(handles.currentYlabel);
title(filename,'Interpreter','none')

% Determine the legend

legend('X','Y')



% Plot fitted steps:
if length(handles.stepVector) == length(handles.currentPlotT)
    if dispLong == 1
        handles.stepLineHandle = plot( ...
            handles.stepVector, '-', 'LineWidth', 2, 'Color', [0.2 0.9 0]);
    end
    if dispShort == 1
        handles.offAxisStepHandle = plot( ...
            handles.shortStepVector, '-', 'LineWidth', 2, 'Color', [0.3 0.3 0.9]);
    end
    if handles.FilteredFlag == 1 && FilteredPlotOffset ~= 0
        if dispLongFilt == 1
            handles.stepLineHandleFilt = plot( ...
                handles.stepVector+FilteredPlotOffset, '-', 'LineWidth', 2, 'Color', [0.3 0.3 0.9]);
        end
        if dispShortFilt == 1
            handles.offAxisStepHandleFilt = plot( ...
                handles.shortStepVector+FilteredPlotOffset, '-', 'LineWidth', 2, 'Color', [0.3 0.3 0.9]);
    
        end
    end
end

% JS Edit 2022/09/27 add in neighbors from other channel
if ~isempty(handles.neighbors)
    % colorlist = summer(round(1.2*length(handles.neighbors)+1)); %shift so that only green
    colorlist = gray(2*length(handles.neighbors)); % only the first half
    colorlist_offaxis = copper(2*length(handles.neighbors)); colorlist_offaxis = flipud(colorlist_offaxis);
    for n = 1:length(handles.neighbors)
        nb = handles.neighbors{n};
        % plot(nb(:,1), nb(:,2), 'Color', colorlist(n,:), 'tag', 'neighbors')
        % plot(nb(:,1), nb(:,3), '--', 'Color', colorlist(n,:), 'tag', 'neighbors')
        scatter(nb(:,1), nb(:,2), 10, 'filled', 'MarkerEdgeColor', [0,0,0], 'MarkerFaceColor', colorlist(n,:), 'tag', 'neighbors')
        scatter(nb(:,1), nb(:,3), 10, 'filled', 'MarkerEdgeColor', [0,0,0], 'MarkerFaceColor', colorlist_offaxis(n,:), 'tag', 'neighbors')
    end
end

% Plot fitted line:
% if length(handles.lineVector) == length(handles.currentPlotT)
%     handles.linfitLineHandle = plot(handles.currentPlotT, ...
%         handles.lineVector, 'g-', 'LineWidth', 2);
% end

hold off


if keepLimits % keep old X and Y limits   
    set(handles.axes1,'Ylim', newYlim);
    set(handles.axes1,'Xlim', newXlim);
end



% Update handles structure
guidata(hObject, handles);