%% 
function varargout = FIONAviewer(varargin)
% m-file file for FIONAviewer.fig
%
%   Written by Vladislav Belyy
%   Last updated on 03/08/2012 
%
% See also: GUIDE, GUIDATA, GUIHANDLES


% Last Modified by GUIDE v2.5 14-Jan-2013 14:47:04
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @FIONAviewer_OpeningFcn, ...
    'gui_OutputFcn',  @FIONAviewer_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before FIONAviewer is made visible.
function FIONAviewer_OpeningFcn(hObject, eventdata, handles, varargin) %#ok
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FIONAviewer (see VARARGIN)

setappdata(0,'UseNativeSystemDialogs',false)

handles.figureHandle = hObject;
% MainPath=pwd;
% SubFunctions=strcat(MainPath,'\SubFns');
% addpath(SubFunctions)
handles.BFlogic=1;



set(handles.FilterData,'string','Filter/Decimate Data?')

handles.display = struct('filtered', 0, 'trapPos', 0, ...
    'gridLines', 0, 'center', 0 ... % 0 = No centering, 1 = median, 2 = user values
    ); % Store displayed elements


handles.rawPSD1Data_X = 0; % unprocessed data from the PSD
handles.rawPSD1Data_Y = 0;
handles.t=0; % time axis
handles.trackAngle = 0; % track angle

% Bead position; formerly handles.X and handles.Y
handles.PSD1Data_Long = 0; % Rotated bead position along track
handles.PSD1Data_Short = 0; % Rotatedposition across track

handles.PSD1Data_Long_Filt = 0; % Filtered bead position
handles.PSD1Data_Short_Filt = 0;
handles.t_Filt = 0; % Filtered/decimated time

% Trap position; formerly handles.Flong and handles.Fshort
handles.trapPosLong=0;
handles.trapPosShort=0;

handles.trapPosLong_Filt=0; % Filtered/decimated trap position
handles.trapPosShort_Filt=0;

                                    
handles.currentPlotT = 0; % Current plot's time axis
handles.currentPlotPSD_Long = 0; % Currently plotted PSD signal, long axis
handles.currentPlotPSD_Short = 0; % Currently plotted PSD signal, short
handles.currentPlotTrap_Long = 0; % Currently plotted trap position, long
handles.currentPlotTrap_Short = 0; % Currently plotted trap position, short
handles.currentXlabel='frame';
handles.currentYlabel='position (nm)'; 

handles.KX=0.05;    % Spring constant, x
handles.KY=0.05;    % Spring constant, y

handles.stepVector = 0; % will store step results
handles.lineVector = 0; % will store line fit results
handles.lineFitCoeffs = []; % Will store line fit coefficients


% Styling 
handles.currentstyleX='b-'; % Currently selected plot style
handles.currentstyleY='r-';
handles.currentstyleTrap='m-';
handles.FilteredFlag = 0;


handles.currentPath='C:\';


% Initialize primary axes
set(handles.axes1,'Ygrid','on')
axes(handles.axes1); %#ok<MAXES>
plot(0,0)
hold on 
xlabel(handles.currentXlabel);
ylabel(handles.currentYlabel);
title('No data selected')
hold off

% % Choose default command line output for FIONAviewer
% handles.output = hObject;
% % Update handles structure
% guidata(hObject, handles);

% JS Edit to make it automatically load what was just run
handles.xydisplayed = 0; %use to keep track whether we are an xy or yx display
if ~isempty(varargin)
    [pathname, filename, ext] = fileparts(varargin(1));
    % run what LoadFile_Callback usually does
    handles = LoadNewDataFile(hObject, handles, fullfile(pathname,'/'), strcat(filename,ext));
    
    keepLimits = 0; % reset Y-limits
    handles = PlotData(hObject, handles, keepLimits); % Plot the data
    
    set(handles.StepsFilename, 'String', fullfile(pathname, strcat(filename,'.mat')));
end

% JS Edit for user preferences related to gridlines, user can change
% if desired
handles.GridDiv.String = '32 nm';
handles.display.gridLines = 1;
set(handles.GridLines,'string','Remove GridLines');
% End of JS Edit

% Choose default command line output for FIONAviewer
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);


% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(~, eventdata, handles)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed


key = eventdata.Key; % determine which key was pressed

% Determine if the user is trying to add or remove steps manually:
manualStepsButton = get(handles.ManualStepFitting, 'String');
manuallyAdjust = strcmp(manualStepsButton, 'Exit');

if manuallyAdjust
       
    if strcmp(key, 't') % toggle between add/delete
        currentString = get(handles.AddDeleteStep, 'String');

        if strcmp(currentString, 'Add')
            set(handles.AddDeleteStep, 'String', 'Delete');
        else
            set(handles.AddDeleteStep, 'String', 'Add');
        end
    end
end






% --- Executes on button press in LoadFile.
function LoadFile_Callback(hObject, eventdata, handles) %#ok



% Ask user to choose file to display
pathcat = get(handles.FilePath,'string');
nameTemplate = [pathcat, '*_fiona.mat'];
[filename, path, filterindex] = uigetfile(nameTemplate, ...
    'Select File to Display (_fiona.mat)'); %#ok

if filename ~= 0 % user selected a name and didn't hit cancel
    handles.xydisplayed = 0;
    handles = LoadNewDataFile(hObject, handles, path, filename);
    
    keepLimits = 0; % reset Y-limits
    handles = PlotData(hObject, handles, keepLimits); % Plot the data
    
end
    
% Update handles structure
guidata(hObject, handles);



function FrameLength_Callback(hObject, ~, handles)

FrameTime = str2double(get(handles.FrameLength, 'string'))/1000; % in sec
handles.t = 0:FrameTime:(length(handles.PSD1Data_Long)-1)*FrameTime;

keepLimits = 0; % reset Y-limits
handles = PlotData(hObject, handles, keepLimits); % Plot the data
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in RePlot.
function RePlot_Callback(hObject, eventdata, handles) %#ok

keepLimits = 0; % Reset limits
handles = PlotData(hObject, handles, keepLimits); % Re-plot data
% Save the handles structure.
guidata(hObject,handles)

% --- Executes on selection change in LinesOrDots. 
function LinesOrDots_Callback(hObject, eventdata, handles) %#ok

keepLimits = 2; % keep previous limits strictly
handles = PlotData(hObject, handles, keepLimits); %Plot data
% Save the handles structure.
guidata(hObject,handles)

% --- Executes on selection change in YScaleMenu.
function YScaleMenu_Callback(hObject, eventdata, handles) %#ok


contents = get(hObject,'String');
currSelection = contents{get(hObject,'Value')};

switch currSelection
    case 'nm and rotated'
        set(handles.GridDiv, 'String', '16 nm');
    case 'Force (pN)'
        set(handles.GridDiv, 'String', '0.5 pN');
    otherwise
        set(handles.GridDiv, 'String', '0.01');
end

keepLimits = 0; % do not keep limits
handles = PlotData(hObject, handles, keepLimits); %Plot data
% Save the handles structure.
guidata(hObject,handles)



% --- Executes on button press in newSpringConstButton.
function newSpringConstButton_Callback(hObject, eventdata, handles) %#ok

% Ask the user for new Kx and Ky:
parastr=inputdlg({'Kx' 'Ky'}, ...
    'Obtaining the Trap Stiffness parameters',1,{'0.05' '0.05'});

if length(parastr) == 2 % user successfully provided two parameters
    
    for i=1:2
        
        paras(i)=str2double(parastr(i)); %#ok
        
    end
    
    set(handles.Kx,'string',num2str(paras(1)));
    set(handles.Ky,'string',num2str(paras(2)));
    
    handles.KX = paras(1);
    handles.KY = paras(2);
    
    keepLimits = 1; % keep previous limits
    handles = PlotData(hObject, handles, keepLimits); %Plot data
    
end

% Save the handles structure.
guidata(hObject,handles)


% --- Executes on button press in ZoomOut.
function ZoomOut_Callback(hObject, eventdata, handles) %#ok

axes(handles.axes1) %#ok<MAXES>
pan off
zoom(0.5)

zoom on



% --- Executes on button press in Pan.
function Pan_Callback(hObject, eventdata, handles) %#ok

axes(handles.axes1) %#ok<MAXES>
panOrZoom=get(handles.Pan,'string');

panOrZoom=panOrZoom(1,:);

if strcmp(panOrZoom, 'Pan')

    zoom off
    pan on

    set(handles.Pan,'string','Zoom');

elseif strcmp(panOrZoom, 'Zoom')

    zoom on

    pan off

    set(handles.Pan,'string','Pan')

end


% --- Executes on button press in SaveCurrentView.
function SaveCurrentView_Callback(hObject, eventdata, handles) %#ok
% Creates an export figure containing a copy of the primary axes

ExportFig = figure('visible','off');
newax = copyobj(handles.axes1, ExportFig);
set(newax, 'units', 'normalized', 'position', [0.1 0.1 0.8 0.8]);
set(ExportFig, 'visible', 'on')
% JS Edit 2022/10/04 to allow for automatic saving of figs (less hastle)
[path,name,~] = fileparts(get(handles.StepsFilename,'String'));
if handles.xydisplayed
    savefig(ExportFig, fullfile(path,strcat(name,'.fig')))
else
    name = name(1:end-6);
    savefig(ExportFig, fullfile(path,strcat(name,'_yx_fiona.fig')))
end
close(ExportFig)



% --- Executes on button press in SaveFile.
function SaveFile_Callback(hObject, eventdata, handles) %#ok
% Saves whatever is currently visible on the main plot as a space-delimited
% file, with the first column storing time in seconds and the remaining
% columns storing whatever else is currently plotted in whatever units it
% is plotted. This info is not saved anywhere in the file, so the user
% should be careful to note the units either in the file name or in their
% lab notebook

currFileName = get(handles.FileName, 'String');
fileMask = [handles.currentPath, '*.txt'];

% Add any additional information to the file name:
addDetails = '_TXYdata.txt';

proposedFileName = strrep(currFileName, 'PSDsignals.txt', addDetails);

% Prompt user for save file name and location
[file,path] = uiputfile(fileMask,'Choose location of new data file', ...
    proposedFileName);

% Create the final file name after the user has had the option to modify it
finalFileName = strcat(path, file);

% Get the data currently plotted on axes1:
children = get(handles.axes1, 'Children');
x = get(children, 'Xdata');
y = get(children, 'Ydata');

% Generate data matrix, with the first column being time and the remaining
% columns being whatever else is currently plotted on the main graph:
M = x{1,1}';
lengthM = length(M);
i = 1;
while length(y{i,1}) == lengthM
    M = [M, y{i,1}']; %#ok<AGROW>
    i = i+1;
end

% Write the data matrix to file with spaces as delimiters:
dlmwrite(finalFileName, M, ' ');

disp('New file created successfully!');

% Save the handles structure.
guidata(hObject,handles)



% --- Executes on button press in GridLines.
function GridLines_Callback(hObject, eventdata, handles)%#ok

Option = get(handles.GridLines,'string');

switch Option

    case 'Add Grid Lines'
        
        handles.display.gridLines = 1;
        set(handles.GridLines,'string','Remove GridLines')

    case 'Remove GridLines'
    
        handles.display.gridLines = 0;
        set(handles.GridLines,'string','Add Grid Lines')

end

keepLimits = 2; % Keep the old limits strictly
handles = PlotData(hObject, handles, keepLimits); % Plot new data

% Save the handles structure.
guidata(hObject,handles)


% --- Executes on button press in DisplayStalls.
function DisplayStalls_Callback(hObject, eventdata, handles) 
%{
x=handles.X;
y=handles.Y;
t=handles.t;

Data=[x' y'];

upper=str2double(get(handles.upper,'string'));
lower=str2double(get(handles.lower,'string'));
%
% cd SubFns

StallStats=display_stalls(100,upper,lower,handles.SR,Data);

FilePath=handles.currentPath;

NewFilePath=nameXYfilefromPSDsignals(FilePath,'Stalls',0);
%
% cd ..

% t=StallStats.T;
% x=StallStats.X;

save(NewFilePath,'-STRUCT','StallStats')

%  fit=StepStats.Fit;
%
% axes(handles.stallsaxes)
%
% ind=1;
%
% cumNP=0;
%
% length(t)
% length(x)
% length(fit)
% StallStats.NumberStalls
% StallStats.NPinStall
%
% for i=1:StallStats.NumberStalls
%
%     NP=StallStats.NPinStall;
%
%     cumNP+1
%     cumNP+NP(i)
%
%     Tpause=t(cumNP+1:cumNP+NP(i));
%     Xpause=x(cumNP+1:cumNP+NP(i));
%     fitpause=fit(cumNP+1:cumNP+NP(i));
%
%     cumNP=cumNP+NP(i);
%
%     plot(Tpause,Xpause,'b.',Tpause,fitpause,'r-')
%
%     hold on
%
% end
% xlabel('time (s)')
% ylabel('position (nm)')
%
% hold off
%
% axes(handles.axes4)
% hist(StepStats.Sizes,20)
% xlabel('step sizes (nm)')

%}




% --- Executes on button press in FilterData.
function FilterData_Callback(hObject, eventdata, handles) %#ok

%handles = ResetStepFitter(handles);

%DisplayType=get(handles.FilterData,'string');

%if strcmp(DisplayType,'Filter/Decimate Data') % Filter data
    
% JS Edit 2023/02/10 just to try and play around with a set filter spec
    %handles = FilterData(hObject, handles); %uncomment if want original
    %FilterData code
    handles = StepThresholdRemoval(hObject, handles);
    % End of JS Edit 2023/02/10

    %set(handles.FilterData,'string','Raw Data?')
    handles.display.filtered = 1;
%end
% else % Display raw data
%     
%     set(handles.FilterData,'string','Filter/Decimate Data?')
%     handles.display.filtered = 0;
% end

keepLimits = 2; % Keep the old limits strictly
handles = PlotData(hObject, handles, keepLimits); % Plot new data

% Save the handles structure.
guidata(hObject,handles)

% --- Executes on selection change in Recording_Type.
function Recording_Type_Callback(hObject, eventdata, handles) %#ok
% Hints: contents = get(hObject,'String') returns Recording_Type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Recording_Type

contents = get(hObject,'String');
currentSel = contents{get(hObject,'Value')};

handles.recordingType = currentSel;

keepLimits = 0; % Reset limits
handles = PlotData(hObject, handles, keepLimits); % Plot new data

% Save the handles structure.
guidata(hObject,handles)

% --- Executes on button press in TrapPositionToggle.
function TrapPositionToggle_Callback(hObject, eventdata, handles) %#ok

status=get(handles.TrapPositionToggle,'string');

if strcmp(status,'Display Trap Position')

    set(handles.TrapPositionToggle,'string','Remove Trap Position');
    handles.display.trapPos = 1;

elseif strcmp(status,'Remove Trap Position')

    set(handles.TrapPositionToggle,'string','Display Trap Position')
    handles.display.trapPos = 0;
    
end

keepLimits = 2; % Keep the old limits strictly
handles = PlotData(hObject, handles, keepLimits); % Plot new data

% Save the handles structure.
guidata(hObject,handles)

% --- Executes on button press in Load_Previous.
function Load_Previous_Callback(hObject, eventdata, handles) %#ok
% Loads the previous PSDsignals file in the current directory

path = get(handles.FilePath, 'string');
file = get(handles.FileName, 'string');

% Obtain a list of all PSDsignals.txt files in the current directory
filenameTemplate = [path, '*_fiona.mat'];
dataFilesStruct = dir(filenameTemplate);
% Convert the structure to cells 
fileNames = cell(length(dataFilesStruct), 1);
for i = 1:length(fileNames)
    fileNames(i) = cellstr(dataFilesStruct(i).name);
end

% Find index of the current file
currentIndex = strmatch(file, fileNames);
% If on yx display, keep same index and switch to xy
% If on xy display, then move to the previous if possible
newIndex = currentIndex - handles.xydisplayed;

% Do not allow the array to exceed bounds
if newIndex > 0
    newFilename = char(fileNames(newIndex));
    % Load the new data file
    handles = LoadNewDataFile(hObject, handles, path, newFilename);

    keepLimits = 0; % reset Y-limits
    handles = PlotData(hObject, handles, keepLimits); % Plot the data
    
    % JS Edit to set for next trace, which is usually a _xy file so save
    % slightly differently
    set(handles.StepsFilename, 'String', fullfile(path, strcat(newFilename)));
%     if handles.xydisplayed
%         set(handles.StepsFilename, 'String', fullfile(path, strcat(newFilename(1:end-10),'_xy.mat')));
%     else
%         set(handles.StepsFilename, 'String', fullfile(path, strcat(newFilename(1:end-10),'_yx.mat')));
%     end
else
    disp('You''ve reached the beginning of the directory');
end
       
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in OpenNext.
function OpenNext_Callback(hObject, eventdata, handles) %#ok

path = get(handles.FilePath, 'string');
file = get(handles.FileName, 'string');

% Obtain a list of all PSDsignals.txt files in the current directory
filenameTemplate = [path, '*_fiona.mat'];
dataFilesStruct = dir(filenameTemplate);
% Convert the structure to cells 
fileNames = cell(length(dataFilesStruct), 1);
for i = 1:length(fileNames)
    fileNames(i) = cellstr(dataFilesStruct(i).name);
end

% Find index of the current file
currentIndex = strmatch(file, fileNames);
% If on xy display, keep same index and switch to yx
% If on yx display, then move to the next if possible
newIndex = currentIndex + ~handles.xydisplayed;

% Do not allow the array to exceed bounds
if newIndex <= length(fileNames)
    newFilename = char(fileNames(newIndex));
    % Load the new data file
    handles = LoadNewDataFile(hObject, handles, path, newFilename);

    keepLimits = 0; % reset Y-limits
    handles = PlotData(hObject, handles, keepLimits); % Plot the data
    % JS Edit to set for next trace, which is usually a _yx file so save
    % slightly differently
    set(handles.StepsFilename, 'String', fullfile(path, strcat(newFilename)));
%     if handles.xydisplayed
%         set(handles.StepsFilename, 'String', fullfile(path, strcat(newFilename(1:end-10),'_xy.mat')));
%     else
%         set(handles.StepsFilename, 'String', fullfile(path, strcat(newFilename(1:end-10),'_yx.mat')));
%     end
else
    disp('You''ve reached the end of the directory');
end
       
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in ShowX.
function ShowX_Callback(hObject, eventdata, handles) %#ok

keepLimits = 1; % keep Y-limits
handles = PlotData(hObject, handles, keepLimits); % Re-plot data
% Save the handles structure.
guidata(hObject,handles)


% --- Executes on button press in ShowY.
function ShowY_Callback(hObject, eventdata, handles) %#ok

keepLimits = 2; % keep Y-limits strictly
handles = PlotData(hObject, handles, keepLimits); % Re-plot data
% Save the handles structure.
guidata(hObject,handles)

function GridDiv_Callback(hObject, eventdata, handles) %#ok
% Hints: get(hObject,'String') returns contents of GridDiv as text
%        str2double(get(hObject,'String')) returns contents of GridDiv as a double

keepLimits = 2; % keep Y-limits strictly
handles = PlotData(hObject, handles, keepLimits); % Re-plot data
% Save the handles structure.
guidata(hObject,handles)


% --- Executes on button press in PrintFilenames.
function PrintFilenames_Callback(~, ~, handles)

path = get(handles.FilePath, 'string');
fileNames = ListFileNames(path);
disp(fileNames);



% --- Executes on button press in PrintPSDCoeffs.
function PrintPSDCoeffs_Callback(~, ~, handles)

path = get(handles.FilePath, 'string');
PrintPSDResponses(path);


% ----- Centering-related callbacks ---------------------------------------
function uipanel9_SelectionChangeFcn(hObject, eventdata, handles) %#ok

% determine new centering parameters
handles = DetermineCentering (hObject, handles);

% Reset step fitter
handles = ResetStepFitter(handles);

keepLimits = 0; % reset Y-limits
handles = PlotData(hObject, handles, keepLimits); % Re-plot data

% Save the handles structure.
guidata(hObject,handles)

function XOffset_Callback(hObject, eventdata, handles) %#ok

% determine new centering parameters
handles = DetermineCentering (hObject, handles);

% Reset step fitter
handles = ResetStepFitter(handles);

keepLimits = 0; % reset Y-limits
handles = PlotData(hObject, handles, keepLimits); % Re-plot data

% Save the handles structure.
guidata(hObject,handles)

function YOffset_Callback(hObject, eventdata, handles) %#ok

% determine new centering parameters
handles = DetermineCentering (hObject, handles);

% Reset step fitter
handles = ResetStepFitter(handles);

keepLimits = 0; % reset Y-limits
handles = PlotData(hObject, handles, keepLimits); % Re-plot data

% Save the handles structure.
guidata(hObject,handles)


%%% Step-fitter callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in FitSteps.
function FitSteps_Callback(hObject, eventdata, handles) %#ok<INUSL>

currSelection = get(hObject, 'String');

if strcmp(currSelection, 'FIT')
    
    %set(hObject, 'String', 'Erase');
    set(handles.ManualStepFitting, 'Visible', 'on');
    set(handles.SaveSteps, 'Visible', 'on');
    set(handles.text72, 'Visible', 'on');
    set(handles.AddDeleteStep, 'Visible', 'on');
    
    % Obtain all the parameters:
    %most of these won't do anything going forward, because SIC don't give
    %a fux
    MaxNumSteps = str2double(get(handles.MaxNumSteps, 'String'));
    ExpSqNoise = str2double(get(handles.ExpSqNoise, 'String'));
    MinDeltaQ = str2double(get(handles.MinDeltaQ, 'String'));
    MinPtsInStep = str2double(get(handles.MinPtsInStep, 'String'));
    
    Xlimits = get(handles.axes1, 'XLim');
    
    %determine start and end points
    startTime =  find(handles.currentPlotT > Xlimits(1), 1);
    endTime = find(handles.currentPlotT > (Xlimits(2)), 1)-1;
    
    if isempty(startTime) 
        startTime = 1; 
    end
    if isempty(endTime)
        endTime = length(handles.currentPlotT);
    end
    
    % make array to fit steps to:
    x = handles.currentPlotPSD_Long(startTime:endTime);
    y = handles.currentPlotPSD_Short(startTime:endTime);
    usage = zeros(1,length(startTime:endTime)) + 1;
    %Note: NaNs removed in parser script
    %x = RemoveNaNs(x);
    
    % fit the steps using my wrapper script to the SIC fitter
    % this version of the wrapper won't fit the off-axis at all, and breaks
    % the trace into subtraces according to usage. See pars_and_fit for
    % more details.
    fitted_trace_6col = parse_and_fit_twosides(x,y,usage,1,500000);
    
    % Build step vectors
    
    %handles.stepVector = NaN(1,length(handles.currentPlotPSD_Long));
    handles.stepVector(startTime:endTime) = (fitted_trace_6col(:,3))';
    handles.shortStepVector(startTime:endTime) = (fitted_trace_6col(:,4))';

else
   handles = ResetStepFitter(handles);
end


keepLimits = 2; % keep Y-limits strictly the same
handles = PlotData(hObject, handles, keepLimits); % Re-plot data

% Save the handles structure.
guidata(hObject,handles)



% This callback gets triggered if user clicks somewhere on the figure and
% is used to add or remove steps manually
function figure1_WindowButtonDownFcn(hObject, ~, handles)

% Determine if the user is trying to add or remove steps manually:
manualStepsButton = get(handles.ManualStepFitting, 'String');
manuallyAdjust = strcmp(manualStepsButton, 'Exit'); 

if manuallyAdjust
    
    % First, determine if the user clicked inside the axes:
    clickCoords = get(handles.axes1, 'CurrentPoint'); % in axes units
    clickX = clickCoords(1,1);
    clickY = clickCoords(1,2);
    
    Xlims = get(handles.axes1, 'XLim');
    Ylims = get(handles.axes1, 'YLim');
    inBounds = clickX>Xlims(1) && clickX<Xlims(2) && clickY>Ylims(1) && ...
        clickY<Ylims(2);
    
    if inBounds
        % Add or remove step:        
        handles = AddRmvStepManually(hObject, handles, clickX);
    end
end
% Save the handles structure.
guidata(hObject,handles)


% --- Executes on button press in AddDeleteStep.
function AddDeleteStep_Callback(hObject, ~, handles)

currentString = get(hObject, 'String');

if strcmp(currentString, 'Add')
    set(hObject, 'String', 'Delete');
else
    set(hObject, 'String', 'Add');
end
% Save the handles structure.
guidata(hObject,handles)

% --- Executes on button press in ManualStepFitting.
function ManualStepFitting_Callback(hObject, ~, handles)

currentString = get(hObject, 'String');

if strcmp(currentString, 'Adjust')
    set(hObject, 'String', 'Exit');
    
    zoom off
    pan off
    
else
    set(hObject, 'String', 'Adjust');
    set(handles.Pan,'string','Pan');
    zoom on
end
% Save the handles structure.
guidata(hObject,handles)



% --- Executes on button press in ChangeStepFilename.
function ChangeStepFilename_Callback(~, ~, handles)

% extract the current directory from filename:
currFullFileName = get(handles.StepsFilename, 'String');
slashPositions = strfind(currFullFileName, '\');
lastSlashPosition = slashPositions(end);
currDir = currFullFileName(1:lastSlashPosition);

filterSpec = [currDir, '*.mat'];
DialogTitle = 'Please select a file to save step-fitting results';


[fileName,pathName,filterIndex] = uiputfile(filterSpec,DialogTitle);  %#ok<NASGU>


if ~isequal(fileName, 0) && ~isequal(pathName, 0)
    set(handles.StepsFilename, 'String', [pathName, fileName]);
end


% --- Executes on button press in SaveSteps.
function SaveSteps_Callback(~, ~, handles)

fileName = get(handles.StepsFilename, 'String');
WriteStepsToFile(handles, fileName);


%%% Line-fitter callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in FitLine.
function FitLine_Callback(hObject, ~, handles)

currSelection = get(hObject, 'String');
if strcmp(currSelection, 'FIT')
    
    set(hObject, 'String', 'Erase');
    set(handles.SaveLineFit, 'Visible', 'on');
       
    Xlimits = get(handles.axes1, 'XLim');
    
    %determine start and end points
    startTime =  find(handles.currentPlotT > Xlimits(1), 1);
    endTime = find(handles.currentPlotT >= (Xlimits(2)), 1)-1;
    
    % select to fit line to:
    x = real(RemoveNaNs(handles.currentPlotPSD_Long(startTime:endTime)));
    t = handles.currentPlotT(startTime:endTime);
    
    % fit to straight line
    linFit = polyfit(t, x, 1);
    handles.lineFitCoeffs = linFit;
    
    % Build line vector
    handles.lineVector = NaN(1,length(handles.currentPlotPSD_Long));
    handles.lineVector(startTime:endTime) = t*linFit(1) + linFit(2);

else % Erase previously fitted line
    handles.lineVector = 0;
    set(handles.FitLine, 'String', 'FIT');
    set(handles.SaveLineFit, 'Visible', 'off');
end

keepLimits = 2; % keep Y-limits strictly the same
handles = PlotData(hObject, handles, keepLimits); % Re-plot data

% Save the handles structure.
guidata(hObject,handles)

function lineFitFilename_Callback(~, ~, handles)

% extract the current directory from filename:
currFullFileName = get(handles.lineFitFilenameEdit, 'String');
slashPositions = strfind(currFullFileName, '\');
lastSlashPosition = slashPositions(end);
currDir = currFullFileName(1:lastSlashPosition);

filterSpec = [currDir, '*.mat'];
DialogTitle = 'Please select a file to save line-fitting results';

[fileName,pathName,filterIndex] = uiputfile(filterSpec,DialogTitle);  %#ok<NASGU>

if ~isequal(fileName, 0) && ~isequal(pathName, 0) % check validity
    set(handles.lineFitFilenameEdit, 'String', [pathName, fileName]);
end


% --- Executes on button press in SaveLineFit.
function SaveLineFit_Callback(~, ~, handles)

fileName = get(handles.lineFitFilenameEdit, 'String');
WriteLineFitToFile(handles, fileName);

%%% Unused calbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% In all of the below callbacks:
% hObject    handle to the object (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function lambda_Callback(hObject, eventdata, handles) 
% Hints: get(hObject,'String') returns contents of lambda as text
%        str2double(get(hObject,'String')) returns contents of lambda as a double

% --- Executes during object creation, after setting all properties.
function lambda_CreateFcn(hObject, eventdata, handles) 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles) 

% --- Executes during object creation, after setting all properties.
function Recording_Type_CreateFcn(hObject, eventdata, handles) 
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CutoffFreq_Callback(hObject, eventdata, handles) 
% Hints: get(hObject,'String') returns contents of CutoffFreq as text
%        str2double(get(hObject,'String')) returns contents of CutoffFreq as a double

% --- Executes during object creation, after setting all properties.
function CutoffFreq_CreateFcn(hObject, eventdata, handles) 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in CentreSignal.
function CentreSignal_Callback(hObject, eventdata, handles) 
% Hint: get(hObject,'Value') returns toggle state of CentreSignal

% --- Executes during object creation, after setting all properties.
function XOffset_CreateFcn(hObject, eventdata, handles) 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function YOffset_CreateFcn(hObject, eventdata, handles) 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles) 
% Hint: get(hObject,'Value') returns toggle state of checkbox7

function FeedbackConversionPara_Callback(hObject, eventdata, handles) 
% Hints: get(hObject,'String') returns contents of FeedbackConversionPara as text
%        str2double(get(hObject,'String')) returns contents of FeedbackConversionPara as a double


% --- Executes during object creation, after setting all properties.
function FeedbackConversionPara_CreateFcn(hObject, eventdata, handles) 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ExpStepSize_Callback(hObject, eventdata, handles) 
% Hints: get(hObject,'String') returns contents of ExpStepSize as text
%        str2double(get(hObject,'String')) returns contents of ExpStepSize as a double

% --- Executes during object creation, after setting all properties.
function ExpStepSize_CreateFcn(hObject, eventdata, handles) 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function upper_Callback(hObject, eventdata, handles) 
% Hints: get(hObject,'String') returns contents of upper as text
%        str2double(get(hObject,'String')) returns contents of upper as a double


% --- Executes during object creation, after setting all properties.
function upper_CreateFcn(hObject, eventdata, handles) 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit9_Callback(hObject, eventdata, handles)  
% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles) 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in PSD1.
function PSD1_Callback(hObject, eventdata, handles) 
% Hint: get(hObject,'Value') returns toggle state of PSD1

% --- Executes on button press in PSD2.
function PSD2_Callback(hObject, eventdata, handles) 
% Hint: get(hObject,'Value') returns toggle state of PSD2

% --- Executes on button press in TrapPos.
function TrapPos_Callback(hObject, eventdata, handles) 
% Hint: get(hObject,'Value') returns toggle state of TrapPos

% --- Executes on button press in BFLogic.
function BFLogic_Callback(hObject, eventdata, handles) 
% Hint: get(hObject,'Value') returns toggle state of BFLogic

function WindowLength_Callback(hObject, eventdata, handles) 
% Hints: get(hObject,'String') returns contents of WindowLength as text
%        str2double(get(hObject,'String')) returns contents of WindowLength as a double

% --- Executes during object creation, after setting all properties.
function WindowLength_CreateFcn(hObject, eventdata, handles) 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit11_Callback(hObject, eventdata, handles) 
% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double

% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles) 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function FilterOrder_Callback(hObject, eventdata, handles) 
% Hints: get(hObject,'String') returns contents of FilterOrder as text
%        str2double(get(hObject,'String')) returns contents of FilterOrder as a double

% --- Executes during object creation, after setting all properties.
function FilterOrder_CreateFcn(hObject, eventdata, handles) 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function GridDiv_CreateFcn(hObject, eventdata, handles) 
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on mouse press over stallsaxes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles) 
get(hObject,'xlim')

disp('You clicked on axis');

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over newSpringConstButton.
function newSpringConstButton_ButtonDownFcn(hObject, eventdata, handles) 

% --- Executes on key press over newSpringConstButton with no controls selected.
function newSpringConstButton_KeyPressFcn(hObject, eventdata, handles) 

% --- Executes during object creation, after setting all properties.
function LinesOrDots_CreateFcn(hObject, eventdata, handles) 
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Outputs from this function are returned to the command line.
function varargout = FIONAviewer_OutputFcn(hObject, eventdata, handles)   %#ok<STOUT>
% varargout  cell array for returning output args (see VARARGOUT);

% --- Executes during object creation, after setting all properties.
function YScaleMenu_CreateFcn(hObject, eventdata, handles) 
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function LowF_Callback(hObject, eventdata, handles) 
% Hints: get(hObject,'String') returns contents of LowF as text
%        str2double(get(hObject,'String')) returns contents of LowF as a double

% --- Executes during object creation, after setting all properties.
function LowF_CreateFcn(hObject, eventdata, handles) %#ok
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Get default command line output from handles structure
varargout{1} = handles.output; %#ok

function edit17_Callback(hObject, eventdata, handles) 
% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in MedianFilter.
function MedianFilter_Callback(hObject, eventdata, handles) 
% Hint: get(hObject,'Value') returns toggle state of MedianFilter

% --- Executes during object creation, after setting all properties.
function Decimate_Factor_CreateFcn(hObject, eventdata, handles) %#ok<*DEFNU,*INUSD>
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function MaxNumSteps_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of MaxNumSteps as text
%        str2double(get(hObject,'String')) returns contents of MaxNumSteps as a double


% --- Executes during object creation, after setting all properties.
function MaxNumSteps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxNumSteps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ExpSqNoise_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of ExpSqNoise as text
%        str2double(get(hObject,'String')) returns contents of ExpSqNoise as a double

% --- Executes during object creation, after setting all properties.
function ExpSqNoise_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function MinDeltaQ_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of MinDeltaQ as text
%        str2double(get(hObject,'String')) returns contents of MinDeltaQ as a double

% --- Executes during object creation, after setting all properties.
function MinDeltaQ_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function MinPtsInStep_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of MinPtsInStep as text
%        str2double(get(hObject,'String')) returns contents of MinPtsInStep as a double


% --- Executes during object creation, after setting all properties.
function MinPtsInStep_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function StepsFilename_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of StepsFilename as text
%        str2double(get(hObject,'String')) returns contents of StepsFilename as a double

% --- Executes during object creation, after setting all properties.
function StepsFilename_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function DecimationFactor_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of DecimationFactor as text
%        str2double(get(hObject,'String')) returns contents of DecimationFactor as a double


% --- Executes during object creation, after setting all properties.
function DecimationFactor_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over PrintFilenames.
function PrintFilenames_ButtonDownFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function NewParameters_CreateFcn(hObject, eventdata, handles)

% --- Executes on mouse press over figure background.
function figure1_ButtonDownFcn(hObject, eventdata, handles)
% disp('You clicked on figure')

% --- Executes during object creation, after setting all properties.
function lineFitFilenameEdit_CreateFcn(hObject, eventdata, handles)

function lineFitFilenameEdit_Callback(hObject, eventdata, handles)








% --- Executes during object creation, after setting all properties.
function FrameLength_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ShowYFilt.
function ShowYFilt_Callback(hObject, ~, handles)

keepLimits = 2; % keep Y-limits strictly
handles = PlotData(hObject, handles, keepLimits); % Re-plot data
% Save the handles structure.
guidata(hObject,handles)
% hObject    handle to ShowYFilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ShowYFilt


% --- Executes on button press in ShowXFilt.
function ShowXFilt_Callback(hObject, ~, handles)
keepLimits = 2; % keep Y-limits strictly
handles = PlotData(hObject, handles, keepLimits); % Re-plot data
% Save the handles structure.
guidata(hObject,handles)
% hObject    handle to ShowXFilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ShowXFilt



function FiltOffsetDistance_Callback(hObject, ~, handles)

 % Re-plot data
% hObject    handle to FiltOffsetDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FiltOffsetDistance as text
%        str2double(get(hObject,'String')) returns contents of FiltOffsetDistance as a double


% --- Executes during object creation, after setting all properties.
function FiltOffsetDistance_CreateFcn(hObject, ~, handles)
% hObject    handle to FiltOffsetDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadmatfile.
function loadmatfile_Callback(hObject, ~, handles)

pathcat = get(handles.FilePath,'string');
[filename, path, filterindex] = uigetfile(pathcat);
if filename ~= 0 % user selected a name and didn't hit cancel
    handles = LoadMatDataFile(hObject, handles, path, filename);
    keepLimits = 0; % reset Y-limits
    handles = PlotData(hObject, handles, keepLimits); % Plot the data
end

set(handles.FitSteps, 'String', 'Erase');
set(handles.ManualStepFitting, 'Visible', 'on');
set(handles.SaveSteps, 'Visible', 'on');
set(handles.text72, 'Visible', 'on');
set(handles.AddDeleteStep, 'Visible', 'on');
    
guidata(hObject, handles);
% hObject    handle to loadmatfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function UsageStart_Callback(hObject, ~, handles)
% hObject    handle to UsageStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of UsageStart as text
%        str2double(get(hObject,'String')) returns contents of UsageStart as a double



% --- Executes during object creation, after setting all properties.
function UsageStart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UsageStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function UsageStopField_Callback(hObject, eventdata, handles)
% hObject    handle to UsageStopField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of UsageStopField as text
%        str2double(get(hObject,'String')) returns contents of UsageStopField as a double


% --- Executes during object creation, after setting all properties.
function UsageStopField_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UsageStopField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function UsageValue_Callback(hObject, eventdata, handles)
% hObject    handle to UsageValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of UsageValue as text
%        str2double(get(hObject,'String')) returns contents of UsageValue as a double


% --- Executes during object creation, after setting all properties.
function UsageValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UsageValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SetUsage.
function SetUsage_Callback(hObject, ~, handles)
% hObject    handle to SetUsage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
start = str2double(get(handles.UsageStart,'String'));
finish = str2double(get(handles.UsageStopField,'String'));
if finish > length(handles.currentPlotPSD_Long)
    finish = length(handles.currentPlotPSD_Long);
end
value = str2double(get(handles.UsageValue,'String'));

handles.usageVector(start:finish) = value;
disp (['setting uage from ' num2str(start) 'to '  num2str(finish) ' to ' num2str(value)' ])
guidata(hObject, handles);







% --- Executes on button press in pushbutton33.
function pushbutton33_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in addUsageAndFit.
function addUsageAndFit_Callback(hObject, eventdata, handles)
%runs the stepFitter, but only for specified usage limits
%vs using the FIT button which fits the entire windowed trace

%first set the usage over the specified region
startTime = str2double(get(handles.UsageStart,'String'));
endTime = str2double(get(handles.UsageStopField,'String'));
value = str2double(get(handles.UsageValue,'String'));

handles.usageVector(startTime:endTime) = value;
disp (['setting uage from ' num2str(startTime) 'to '  num2str(endTime) ' to ' num2str(value)' ])

%then fit the trace, this is mostly taken from the "FIT" button callback


    set(handles.ManualStepFitting, 'Visible', 'on');
    set(handles.SaveSteps, 'Visible', 'on');
    set(handles.text72, 'Visible', 'on');
    set(handles.AddDeleteStep, 'Visible', 'on');
    
    
    % make array to fit steps to:
    x = handles.currentPlotPSD_Long(startTime:endTime);
    y = handles.currentPlotPSD_Short(startTime:endTime);
    usage = handles.usageVector(startTime:endTime);
    %Note: NaNs removed in parser script
    %x = RemoveNaNs(x);
    
    % fit the steps using my wrapper script to the SIC fitter
    % this version of the wrapper won't fit the off-axis at all, and breaks
    % the trace into subtraces according to usage. See parse_and_fit for
    % more details.
    
    fitted_trace_6col = parse_and_fit_twosides(x,y,usage,1,500000);
    
    % Build step vectors
    
    %handles.stepVector = NaN(1,length(handles.currentPlotPSD_Long));
    handles.stepVector(startTime:endTime) = (fitted_trace_6col(:,3))';
    handles.shortStepVector(startTime:endTime) = (fitted_trace_6col(:,4))';
    keepLimits = 2; % keep Y-limits strictly the same
    handles = PlotData(hObject, handles, keepLimits); % Re-plot data
    
guidata(hObject, handles);


% hObject    handle to addUsageAndFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on SetUsage and none of its controls.

% hObject    handle to SetUsage (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton35.
% this is the "fit only" button
function pushbutton35_Callback(hObject, ~, handles)
%runs the stepFitter, but only for specified usage limits
%vs using the FIT button which fits the entire windowed trace

%first set the usage over the specified region
startTime = str2double(get(handles.UsageStart,'String'));
endTime = str2double(get(handles.UsageStopField,'String'));

disp (['fitting trace between ' num2str(startTime) 'and '  num2str(endTime) ])

%then fit the trace, this is mostly taken from the "FIT" button callback


    set(handles.ManualStepFitting, 'Visible', 'on');
    set(handles.SaveSteps, 'Visible', 'on');
    set(handles.text72, 'Visible', 'on');
    set(handles.AddDeleteStep, 'Visible', 'on');
  
    
    % make array to fit steps to:
    x = handles.currentPlotPSD_Long(startTime:endTime);
    y = handles.currentPlotPSD_Short(startTime:endTime);
    usage = ones(1,length(x));
    %Note: NaNs removed in parser script
    %x = RemoveNaNs(x);
    
    % fit the steps using my wrapper script to the SIC fitter
    % this version of the wrapper won't fit the off-axis at all, and breaks
    % the trace into subtraces according to usage. See parse_and_fit for
    % more details.
    
    fitted_trace_6col = parse_and_fit_twosides(x,y,usage,1,500000);
    
    % Build step vectors
    
    %handles.stepVector = NaN(1,length(handles.currentPlotPSD_Long));
    handles.stepVector(startTime:endTime) = (fitted_trace_6col(:,3))';
    handles.shortStepVector(startTime:endTime) = (fitted_trace_6col(:,4))';
    keepLimits = 2; % keep Y-limits strictly the same
    handles = PlotData(hObject, handles, keepLimits); % Re-plot data
    
guidata(hObject, handles);

% hObject    handle to pushbutton35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
