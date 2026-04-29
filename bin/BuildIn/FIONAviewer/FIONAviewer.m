function app = FIONAviewer(varargin)
    % Modern UIFIGURE-based recreation of the original FIONAviewer GUI
    % Works in MATLAB R2020+ (tested R2023, R2024)
    
    % make sure to add dependencies if launching through Fiesta
    addpath('bin/BuildIn/FIONAviewer/SubFns/')
    % otherwise do a simpler version if already navigated to parent
    % FIONAviewer
    % addpath('SubFns/')

    app.fig = uifigure('Name','FIONAviewer', ...
                       'Position',[100 100 1400 800]);
    
    %% Initialize all necessary fields (in their subfields if possible)
    % We will use subcategorization
    % app.LeftPanel
    app.LeftPanel = [];
        app.LeftPanelFields.StepsFilename = [];
        app.LeftPanelFields.BtnLoadMat = [];
        app.LeftPanelFields.BtnLoadFile = [];
        app.LeftPanelFields.BtnPrev = [];
        app.LeftPanelFields.BtnNext = [];
        app.LeftPanelFields.BtnSave = [];
    
    app.FitPanel = [];
    % app.EditNumSteps = [];
    % app.EditNoise = [];
    % app.EditMinDelta = [];
        app.FitPanelFields.DropFitAxis = [];
        app.FitPanelFields.DropFitString = [];
        app.FitPanelFields.BtnFit = [];
        app.FitPanelFields.BtnAdd = [];
        app.FitPanelFields.BtnDelete = [];
        app.FitPanelFields.MeanSampleWindow = [];
        app.FitPanelFields.MinDeltaChange = [];
        app.FitPanelFields.BtnAutoAlign = [];
    % app.PanI = panInteraction; %Not functional will try again later
    % app.ZoomI = zoomInteraction;
    
    app.FilterPanel = [];
        app.FilterPanelFields.FilterGroup = [];
        app.FilterPanelFields.BtnFilter = [];
        app.FilterPanelFields.FilterWindow = [];
        app.FilterPanelFields.FilterOrder = [];
        app.FilterPanelFields.FilterLambdaPenalty = [];

    app.DisplayPanel = [];
        app.DisplayPanelFields.DropLineStyle = [];
        app.DisplayPanelFields.BtnRefresh = [];
        app.DisplayPanelFields.BtnZoomOut = [];
        app.DisplayPanelFields.ChkGrid = [];

    % app.StepPanel = [];
    % app.EditStart = [];
    % app.EditEnd = [];
    % app.EditValue = [];
    % app.BtnFitStepPanel2 = []; % renamed to avoid overwriting
    % app.BtnAddFit = [];
    % Any data fields used later
    app.Data.t = [];
    app.Data.PSD1Data_Long = [];
    app.Data.PSD1Data_Short = [];
    app.Data.stepVector = [];
    app.Data.shortstepVector = [];
    app.Data.FilterData = [];
    app.Data.FilteredFlag = 0;
    app.Data.StepChangesTable = [];
    app.time_bool = 0;
    app.display = struct();
    
    % this is the storage space for what we will eventually export in "Save"
    app.data = [];

    app.FileName = [];
    app.FilePath = [];

    %% ---------------------------
    %  AXES (Main Plot Area)
    % ----------------------------
    app.AxMain = uiaxes(app.fig, ...
        'Position',[350 100 1000 650], ...
        'Tag','MainAxes');
    title(app.AxMain,'');
    xlabel(app.AxMain,'frame');
    ylabel(app.AxMain,'position (nm)');
    grid(app.AxMain,'on');
    set(app.AxMain, 'FontName', 'Arial', 'FontSize', 10, 'TickDir', 'out', 'LineWidth', 1, 'Box', 'off', 'XColor', 'k', 'YColor', 'k');

    %% -------------------------------------------
    % LEFT PANEL — File Loading + Status
    % -------------------------------------------
    app.LeftPanel = uipanel(app.fig, ...
        'Title','Current File', ...
        'Position',[30 530 300 220]);

    app.LeftPanelFields.StepsFilename = uilabel(app.LeftPanel, ...
        'Position',[20 170 250 30], ...
        'Text','No file loaded');

    app.LeftPanelFields.BtnLoadMat = uibutton(app.LeftPanel, ...
        'Position',[20 135 250 35], ...
        'Text','Load New .mat File', ...
        'ButtonPushedFcn',@(src,evt)onLoadMat(app, []));

    app.LeftPanelFields.BtnLoadFile = uibutton(app.LeftPanel, ...
        'Position',[20 95 250 35], ...
        'Text','Load New File?', ...
        'ButtonPushedFcn',@(src,evt)onLoadNew(app));

    app.LeftPanelFields.BtnPrev = uibutton(app.LeftPanel, ...
        'Position',[20 50 120 35], ...
        'Text','←  Previous', ...
        'ButtonPushedFcn',@(src,evt)onPrev(app, src));

    app.LeftPanelFields.BtnNext = uibutton(app.LeftPanel, ...
        'Position',[150 50 120 35], ...
        'Text','Next  →', ...
        'ButtonPushedFcn',@(src,evt)onNext(app, src));

    app.LeftPanelFields.BtnExport = uibutton(app.LeftPanel, ...
        'Position',[20 10 120 35], ...
        'Text','Export Figure', ...
        'ButtonPushedFcn', @(src,event) onExportFigure(app));

    app.LeftPanelFields.BtnSave = uibutton(app.LeftPanel, ...
        'Position',[150 10 120 35], ...
        'Text','Save', ...
        'ButtonPushedFcn',@(src,evt)onSave(app));

    %% -------------------------------------------
    % FIT STEPS PANEL (Left-mid)
    % -------------------------------------------
    app.FitPanel = uipanel(app.fig, ...
        'Title','Fit steps', ...
        'Position',[20 390 300 120]);

    % uilabel(app.FitPanel,'Position',[10 120 150 20],'Text','Max # of steps:');
    % app.EditNumSteps = uieditfield(app.FitPanel,'numeric', ...
    %     'Position',[160 120 50 22], ...
    %     'Value',20);
    % 
    % uilabel(app.FitPanel,'Position',[10 90 150 20],'Text','Expected sq noise:');
    % app.EditNoise = uieditfield(app.FitPanel,'numeric', ...
    %     'Position',[160 90 50 22], ...
    %     'Value',2);
    % 
    % uilabel(app.FitPanel,'Position',[10 60 150 20],'Text','Min DeltaQ:');
    % app.EditMinDelta = uieditfield(app.FitPanel,'numeric', ...
    %     'Position',[160 60 50 22], ...
    %     'Value',2);

    uilabel(app.FitPanel,'Position',[10 75 60 20],'Text','Fit Axis:');
    app.FitPanelFields.DropFitAxis = uidropdown(app.FitPanel, ...
        'Position',[10 55 120 22], ...
        'Items',{'Long-axis','Short-axis'}, ...
        'Value','Long-axis');

    % --- Create edit field left of the Fit Axis dropdown
    app.FitPanelFields.DropFitString = uieditfield(app.FitPanel, 'text');
    app.FitPanelFields.DropFitString.Position = [140 55 70 22]; % set proper coordinates
    app.FitPanelFields.DropFitString.Tooltip = 'User string for step changes';
    app.FitPanelFields.DropFitString.Value = 'user';


    app.FitPanelFields.BtnFit = uibutton(app.FitPanel, ...
        'Position',[10 15 120 30], ...
        'Text','Fit', ...
        'ButtonPushedFcn',@(src,evt)onFit(app));
    
    app.FitPanelFields.BtnAdd = uibutton(app.FitPanel, ...
        'Position',[150 15 60 30], ...
        'Text','Add', ...
        'ButtonPushedFcn',@(src,evt)onAdd(app, src));
    app.FitPanelFields.BtnAdd.UserData.pressed = false;
    app.FitPanelFields.BtnAdd.BackgroundColor = [1 1 1];
    
    app.FitPanelFields.BtnDelete = uibutton(app.FitPanel, ...
        'Position',[220 15 60 30], ...
        'Text','Delete', ...
        'ButtonPushedFcn',@(src,evt)onDelete(app, src));
    app.FitPanelFields.BtnDelete.UserData.pressed = false;
    app.FitPanelFields.BtnDelete.BackgroundColor = [1 1 1];


    app.FitPanelFields.BtnAutoAlign = uibutton(app.FitPanel, ...
    'Position',[220 55 70 30], ...
    'Text','Align', ...
    'ButtonPushedFcn',@(src,evt)onAlign(app));
    
    % uilabel(app.FitPanel,'Position',[200 50 40 30],'Text','Window:');
    % app.FitPanelFields.MeanSampleWindow = uieditfield(app.FitPanel,'numeric', ...
    %     'Position',[150 55 50 25],'Value', 10);
    
    % uilabel(app.FitPanel,'Position',[140 75 80 30],'Text','Min Change:');
    % app.FitPanelFields.MinDeltaChange = uieditfield(app.FitPanel,'numeric', ...
    %     'Position',[150 55 60 25],'Value', 4);
    
    
    %% MANUAL CLICKS
    app.AxMain.ButtonDownFcn = @(src, event) AxMainClick(app, src, event);
    % HOTKEYS
    app.fig.WindowKeyPressFcn = @(fig,event) onKeyPress(app, event);
    % app.AxMain.ButtonDownFcn = @(src,event) uifigureFocus(app);

    %% Shortcut accelerators (Define all strings hear for continuity throughout the code)
    % User defined shortcuts here
    % app.shortcuts = dictionary("zoom","z", "zoom_out","x", "pan","v", ...
    %     "add","a", "delete","d", "filter","f", "save","s", "export","c");

    % Add menus with Accelerators
    HKmenu = uimenu('Parent',app.fig,'Label','Hot Keys');
    uimenu('Parent',HKmenu,'Label','Zoom','Accelerator','z','Callback',@(src,evt)zoom(app.fig,'toggle'));
    uimenu('Parent',HKmenu,'Label','Reset','Accelerator','x','Callback',@(src,evt)zoom(app.fig,'out'));
    uimenu('Parent',HKmenu,'Label','Pan','Accelerator','v','Callback',@(src,evt)pan(app.fig,'toggle'));

    % uimenu('Parent',HKmenu,'Label','Save','Accelerator','s','Callback',@(src,evt)onSave(app));
    % uimenu('Parent',HKmenu,'Label','Export Fig','Accelerator','c','Callback',@(src,evt)onExportFigure(app));

    %% -------------------------------------------
    % FILTER PANEL (Bottom-left)
    % -------------------------------------------
    app.FilterPanel = uipanel(app.fig, ...
        'Title','Filter Method', ...
        'Position',[20 50 300 180]);

    app.FilterPanelFields.FilterGroup = uibuttongroup(app.FilterPanel, ...
        'Position',[10 50 280 100]);

    uiradiobutton(app.FilterPanelFields.FilterGroup,'Text','Decimate','Position',[10 80 150 20]);
    uiradiobutton(app.FilterPanelFields.FilterGroup,'Text','Median Filter','Position',[10 60 150 20]);
    uiradiobutton(app.FilterPanelFields.FilterGroup,'Text','Running Mean','Position',[10 40 150 20]);
    uiradiobutton(app.FilterPanelFields.FilterGroup,'Text','Butterworth','Position',[10 20 150 20]);
    uiradiobutton(app.FilterPanelFields.FilterGroup,'Text','L1 Piecewise Constant','Position',[10 0 200 20]);
    
    uilabel(app.FilterPanel,'Position',[145 120 80 25],'Text','Dec/Window:');
    app.FilterPanelFields.FilterWindow = uieditfield(app.FilterPanel,'numeric', ...
        'Position',[225 120 50 25],'Value', 10);
    
    uilabel(app.FilterPanel,'Position',[145 93 80 25],'Text','Filter Order:');
    app.FilterPanelFields.FilterOrder = uieditfield(app.FilterPanel,'numeric', ...
        'Position',[225 93 50 25],'Value', 4);
    
    uilabel(app.FilterPanel,'Position',[145 66 80 25],'Text','Lambda:');
    app.FilterPanelFields.FilterLambdaPenalty = uieditfield(app.FilterPanel,'numeric', ...
        'Position',[225 66 50 25],'Value', 15);

    app.FilterPanelFields.BtnFilter = uibutton(app.FilterPanel, ...
        'Position',[150 10 130 35], ...
        'Text','Filter', ...
        'ButtonPushedFcn',@(src,evt)onFilter(app, src));
    app.FilterPanelFields.BtnFilter.UserData.pressed = false;
    app.FilterPanelFields.BtnFilter.BackgroundColor = [1 1 1];

    %% -------------------------------------------
    % DISPLAY CONTROL PANEL (Bottom middle)
    % -------------------------------------------
    app.DisplayPanel = uipanel(app.fig, ...
        'Title','Display control', ...
        'Position',[20 250 300 120]); %'Position',[350 20 450 130]);

    uilabel(app.DisplayPanel,'Position',[10 70 60 20],'Text','Plot marker:');
    app.DisplayPanelFields.DropLineStyle = uidropdown(app.DisplayPanel, ...
        'Position',[80 70 100 22], ...
        'Items',{'Lines','Dots','Both','None'}, ...
        'Value','Lines');

    app.DisplayPanelFields.BtnRefresh = uibutton(app.DisplayPanel, ...
        'Position',[200 65 80 30], ...
        'Text','Refresh', ...
        'ButtonPushedFcn',@(src,evt)onReplot(app));

    % app.DisplayPanelFields.BtnZoomOut = uibutton(app.DisplayPanel, ...
    %     'Position',[200 25 80 30], ...
    %     'Text','Zoom Out', ...
    %     'ButtonPushedFcn',@(s,e)onZoomOut(app));

    app.DisplayPanelFields.ChkGrid = uicheckbox(app.DisplayPanel, ...
        'Position',[10 25 120 30], ...
        'Text','Grid Spacing','Value', 1);

    % uilabel(app.DisplayPanelFields.GridSpacingText,'Position',[10 120 150 20],'Text','Grid Spacing (nm)');
    app.DisplayPanelFields.GridSpacing = uieditfield(app.DisplayPanel,'numeric', ...
        'Position',[100 25 40 30], ...
        'Value', 16);

    % %% -------------------------------------------
    % % STEP FITTER PANEL (Bottom-right, Yellow)
    % % -------------------------------------------
    % app.StepPanel = uipanel(app.fig, ...
    %     'Title','Set Usage For Step Fitter', ...
    %     'Position',[820 20 350 130]);
    % 
    % uilabel(app.StepPanel,'Position',[10 70 60 20],'Text','Start');
    % uilabel(app.StepPanel,'Position',[10 40 60 20],'Text','End');
    % uilabel(app.StepPanel,'Position',[10 10 60 20],'Text','Value');
    % 
    % app.EditStart = uieditfield(app.StepPanel,'numeric','Position',[70 70 60 22]);
    % app.EditEnd   = uieditfield(app.StepPanel,'numeric','Position',[70 40 60 22]);
    % app.EditValue = uieditfield(app.StepPanel,'numeric','Position',[70 10 60 22]);
    % 
    % app.BtnAdd     = uibutton(app.StepPanel,'Text','ADD', ...
    %     'Position',[150 70 80 25], ...
    %     'ButtonPushedFcn',@(s,e)onAdd(app));
    % app.BtnFit     = uibutton(app.StepPanel,'Text','FIT', ...
    %     'Position',[150 40 80 25], ...
    %     'ButtonPushedFcn',@(s,e)onFit(app));
    % app.BtnAddFit  = uibutton(app.StepPanel,'Text','ADD AND FIT', ...
    %     'Position',[150 10 120 25], ...
    %     'ButtonPushedFcn',@(s,e)onAddAndFit(app));
    

    if ~isempty(varargin)
        setappdata(app.fig, 'app', app);
        app = onLoadMat(app, varargin{1});
    end

    % at the end save it all
    setappdata(app.fig, 'app', app);

end

% Load existing MATLAB file
function app = onLoadMat(app, name)
    app = getappdata(app.fig, 'app');
    if isempty(name)
        [file, path] = uigetfile('*.mat');
    else
        [path, file, ~] = fileparts(name);
    end

    if isequal(file,0)
        return;
    end

    % Call your loader
    app = LoadMatDataFile(app, path, file);

    % Update plot on the main axes
    app = updateMainPlot(app);
    setappdata(app.fig, 'app', app);
end

% Load new and convert MATLAB file
function app = onLoadNew(app)
    app = getappdata(app.fig, 'app');
    [file, path] = uigetfile('*.mat');

    if isequal(file,0)
        return;
    end

    % Call your newly modernized loader
    app = LoadNewDataFile(app, path, file);

    % Update plot on the main axes
    app = updateMainPlot(app);
    setappdata(app.fig, 'app', app);
end

function onNext(app, source)
    app = getappdata(app.fig, 'app');
    app = ScrollPrevNext(app, source);
    setappdata(app.fig,'app',app);
end

function onPrev(app, source)
    app = getappdata(app.fig, 'app');
    app = ScrollPrevNext(app, source);
    setappdata(app.fig,'app',app);
end


function onExportFigure(app)
    ExportFigure(app);
end


function app = onSave(app)
    app = getappdata(app.fig, 'app');
    
    data = app.data;
    %this will now re-write traces so be careful!
    % app.LeftPanelFields.StepsFilename.Text
    disp([ 'writing current trace to ' app.LeftPanelFields.StepsFilename.Text ])
    save(app.LeftPanelFields.StepsFilename.Text ,'data');


    setappdata(app.fig, 'app', app);
end


function app = onFit(app)
    app = getappdata(app.fig, 'app');
    if strcmp(app.FitPanelFields.DropFitAxis.Value,'Long-axis')
        lineObj = findobj(app.AxMain, 'Type', 'Line', 'Tag', 'hLongAxis');
        stepFit = SICstepFinder(lineObj.YData);
        % reformat into size N so easily writeable later
        if length(stepFit.StepFit) ~= length(app.Data.t)
            out = repelem(stepFit.StepFit, app.FilterPanelFields.FilterWindow.Value);
            stepFit.StepFit = out(1:length(app.Data.t));
        end
        app.Data.stepVector = stepFit.StepFit; % Update step vector
        cp = find(abs(app.Data.stepVector(2:end) - app.Data.stepVector(1:end-1)));
        app = UserStepChangesTable(app, cp+1, 'fit');
    elseif strcmp(app.FitPanelFields.DropFitAxis.Value,'Short-axis')
        lineObj = findobj(app.AxMain, 'Type', 'Line', 'Tag', 'hShortAxis');
        stepFit = SICstepFinder(lineObj.YData);
        % reformat into size N so easily writeable later
        if length(stepFit.StepFit) ~= length(app.Data.t)
            out = repelem(stepFit.StepFit, app.FilterPanelFields.FilterWindow.Value);
            stepFit.StepFit = out(1:length(app.Data.t));
        end
        app.Data.shortstepVector = stepFit.StepFit; % Update step vector
        cp = find(abs(app.Data.shortstepVector(2:end) - app.Data.shortstepVector(1:end-1)));
        app = UserStepChangesTable(app, cp+1, 'fit');
    end
    app = updateMainPlot(app); %and replot
    % Then let buttons be visible again
    set(app.FitPanelFields.BtnAdd, 'Visible', 'on');
    set(app.FitPanelFields.BtnDelete, 'Visible', 'on');
    set(app.LeftPanelFields.BtnSave, 'Visible', 'on');
    % Get trace ready for export
    app = PackageTrace(app, stepFit);
    % app.data.trace   % Can check that it worked here
    % app.data.trace_yx
    setappdata(app.fig, 'app', app);
end


function app = onAdd(app, btn)
    app = getappdata(app.fig, 'app');
    
    % release other button (hack so that it will just toggle)
    app.FitPanelFields.BtnDelete.BackgroundColor = [1 1 1];
    app.FitPanelFields.BtnDelete.UserData.pressed = false;
    
    % Now decide if you want to toggle this one
    if ~isfield(btn.UserData, 'pressed') || ~btn.UserData.pressed
        % Button is now pressed
        btn.BackgroundColor = [0.92 1 1];  % light cyan
        btn.UserData.pressed = true;
    else
        % Button is now released
        btn.BackgroundColor = [1 1 1];      % default
        btn.UserData.pressed = false;
    end
    setappdata(app.fig, 'app', app);
end


function app = onDelete(app, btn)
    app = getappdata(app.fig, 'app');

    % release other button (hack so that it will just toggle)
    app.FitPanelFields.BtnAdd.BackgroundColor = [1 1 1];
    app.FitPanelFields.BtnAdd.UserData.pressed = false;

    if ~isfield(btn.UserData, 'pressed') || ~btn.UserData.pressed
        % Button is now pressed
        btn.BackgroundColor = [0.92 1 1];  % light cyan
        btn.UserData.pressed = true;
    else
        % Button is now released 
        btn.BackgroundColor = [1 1 1];      % default
        btn.UserData.pressed = false;
    end
    setappdata(app.fig, 'app', app);
end


function app = onAlign(app)
    app = getappdata(app.fig, 'app');
    
    app = AlignToMaxDx(app);
    app = updateMainPlot(app);
    
    setappdata(app.fig, 'app', app);
end



function app = onReplot(app)
    app = getappdata(app.fig, 'app');
    app = updateMainPlot(app);
    setappdata(app.fig, 'app', app);
end


function app = onFilter(app, btn)
    app = getappdata(app.fig, 'app');
    
    if ~isfield(btn.UserData, 'pressed') || ~btn.UserData.pressed
        % Button is now pressed
        btn.BackgroundColor = [0.92 1 1];  % light cyan
        btn.UserData.pressed = true;
        app = FilterData(app);
    else
        % Button is now released 
        btn.BackgroundColor = [1 1 1];      % default
        btn.UserData.pressed = false;
        app.Data.FilteredFlag = 0;
    end
    
    app = updateMainPlot(app);
    % app = updateMainPlot(app);
    setappdata(app.fig, 'app', app);
end


% Main click callback function
function [nearestX, nearestY] = AxMainClick(app, src, event)
    % disp('click fired')
    app = getappdata(app.fig, 'app');
    figure(app.fig) % 1. Focus the figure
    
        % 2. Then find nearest point
    % Check if zoom OR pan mode is active
    if strcmp(zoom(app.fig).Enable, 'on') || strcmp(pan(app.fig).Enable, 'on')
        return   % <-- Skip the click behavior
    end
    
    % Otherwise run your nearest-point logic
    cp = event.IntersectionPoint;     % [x y z]
    xclick = cp(1);
    yclick = cp(2);
    
    ax = app.AxMain;
    xRange = diff(ax.XLim);
    yRange = diff(ax.YLim);
    
    % Find the line (if only one line, just pick it)
    if strcmp(app.FitPanelFields.DropFitAxis.Value,'Long-axis')
        lineObj = findobj(ax, 'Type', 'Line', 'Tag', 'hFitLine'); %line to find nearest point to Long-axis
    else
        lineObj = findobj(ax, 'Type', 'Line', 'Tag', 'hFitShortLine'); %line to find nearest point to Short-axis fit
    end

    % If there is nothing there, just leave the function
    if isempty(lineObj)
        return
    end
    
    xline = lineObj.XData;
    yline = lineObj.YData;
    
    % Normalize both click and line points
    xdataN = (xline - ax.XLim(1)) / xRange;
    ydataN = (yline - ax.YLim(1)) / yRange;
    xclickN = (xclick - ax.XLim(1)) / xRange;
    yclickN = (yclick - ax.YLim(1)) / yRange;
    
    % Vectorized distance computation
    dist = hypot(xdataN - xclickN, ydataN - yclickN);
    
    [~, idx] = min(dist);
    
    nearestX = xline(idx);
    nearestY = yline(idx);

    % FOR DEBUGGING
    % Display the selected point (example action)
    % fprintf('Nearest point: (%.3f, %.3f)\n', nearestX, nearestY);

    % Optional: mark it 
    % hold(app.AxMain, 'on');
    % if isfield(app.DisplayPanelFields,'ClickMarker') && isvalid(app.DisplayPanelFields.ClickMarker)
    %     delete(app.DisplayPanelFields.ClickMarker)
    % end
    % app.DisplayPanelFields.ClickMarker = plot(app.AxMain, nearestX, nearestY, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
    % hold(app.AxMain, 'off');
    
    % run a function (add/rm steps)
    if app.FitPanelFields.BtnAdd.UserData.pressed || app.FitPanelFields.BtnDelete.UserData.pressed
        app = AddRmvStepManually(app, nearestX, nearestY);
    end
    app = updateMainPlot(app);
    setappdata(app.fig, 'app', app);
end


% HOTKEYS
% function uifigureFocus(app)
%     figure(app.fig); % give figure focus so that it knows to change something
% end
function onKeyPress(app, event)
    key = event.Key;  % string of the key pressed, e.g., 'a', 'd', 'leftarrow'
    mods = event.Modifier; % cell array of modifiers, e.g., {'control'}
    
    % switch between control for pc and command for mac
    if ispc && ismember('control', mods) || ismac && ismember('command', mods)
        switch key
            case 's'
                % disp('Ctrl+s pressed')
                onSave(app)

            case 'c'
                % disp('Ctrl+e pressed')
                onExportFigure(app)

            % case 'v'
            %     % disp('Ctrl+v pressed')
            %     pan(app.fig,'toggle')
            % 
            % case 'z'
            %     % disp('Ctrl+z pressed')
            %     zoom(app.fig,'toggle')
            % 
            % case 'x'
            %     % disp('Ctrl+x pressed')
            %     zoom(app.fig,'out')

        end

    else
        switch key
            case 'a' %a trigger onAdd button
                onAdd(app, app.FitPanelFields.BtnAdd);
            
            case 'd' %d trigger onDelete button
                onDelete(app, app.FitPanelFields.BtnDelete);

            % case 'f' %f trigger onFilter button
            %     onFilter(app, app.FilterPanelFields.BtnFilter);

            % case 'v' %v triggers pan
            %     pan(app.fig,'toggle') %toggle only will turn on, we can't access this in the active figure mode state
            % 
            % case 'z' %z triggers zoom in
            %     zoom(app.fig,'toggle')
            % 
            % case 'x' %x triggers reset (zoom out)
            %     zoom(app.fig,'out')

            % You can add more hotkeys here
        end
    end
end

