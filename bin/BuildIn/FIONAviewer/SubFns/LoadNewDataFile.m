function handles = LoadNewDataFile(hObject, handles, path, filename)
% Reads data from a new file given by filename at the given path and saves 
% all the appropriate new variables to the handles structure
%
% Written by Vladislav Belyy
% Last updated on 11/18/2011


%% Prepare for data reading

tic % debug purposes only - optimizing load time

% Turn off filtering
set(handles.FilterData,'string','Filter/Decimate Data?')
handles.display.filtered = 0;

handles = ResetStepFitter(handles);

% Reset line fitter
% handles.lineVector = 0;
% set(handles.FitLine, 'String', 'FIT');
% set(handles.SaveLineFit, 'Visible', 'off');

DataFilePath=strcat(path,filename);
handles.currentPath=DataFilePath;





%% Read data

Data=ReadDataFile(DataFilePath); % fast data read

tElapsed = toc; % debug purposes only - optimizing time


set(handles.FileName,'string',filename)
set(handles.FilePath,'string',path)



%% Load data 

%creates empty data for four of the columns in my 6-column trace files

% Rotate and convert raw PSD1 signals: creates handles.PSD1Data_Long,
% PSD1Data_Short, and handles.t
handles.PSD1Data_Long = Data(:,1)';
handles.PSD1Data_Short = Data(:,2)';
limit = length(Data(:,1));
handles.stepVector(1:limit) = NaN;
handles.shortStepVector(1:limit) = NaN;
% default usage instructions are zero for data loaded from two-column text
handles.usageVector = 0;
handles.usageVector(1,1:length(Data(:,1))) = 0;
%also initialize a changepoints vector, for consistency
handles.changepoints = 0;
handles.changepoints(1,1:length(Data(:,1))) = 0;

%initialize time vector as just the indices (if its even still used)
handles.t = 1:length(Data(:,1));


% debug purposes only:
disp(['Loaded in: ', num2str(tElapsed), ' seconds']);

% Update handles structure
guidata(hObject, handles);