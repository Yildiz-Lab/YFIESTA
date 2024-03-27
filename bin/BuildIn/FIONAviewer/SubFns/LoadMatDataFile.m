function handles = LoadMatDataFile(~, handles, path, filename)

tic % debug purposes only - optimizing load time

% Turn off filtering
set(handles.FilterData,'string','Filter/Decimate Data?');
handles.display.filtered = 0;

handles = ResetStepFitter(handles);

% Reset line fitter
handles.lineVector = 0;
%set(handles.FitLine, 'String', 'FIT');
%set(handles.SaveLineFit, 'Visible', 'off');

DataFilePath=strcat(path,filename);
handles.currentPath=DataFilePath;

% JS Edit 2024/03/18 for new mat way to save data
trace_curr = load(DataFilePath);
trace_curr = trace_curr.data;

% trace_curr = load (DataFilePath,'-mat','trace');

handles.PSD1Data_Long = trace_curr.trace(:,1)';
handles.PSD1Data_Short = trace_curr.trace(:,2)';
handles.stepVector = trace_curr.trace(:,3)';
handles.shortStepVector = trace_curr.trace(:,4)';
%make a new handle for the usage from the mat file...so we can save it
%right later.s
handles.usageVector = trace_curr.trace(:,6)';

%FrameTime = str2double(get(handles.FrameLength, 'string'))/1000; % in sec
handles.t = 1:length(trace_curr.trace(:,1));

% JS Edit 2024/03/18 for new mat way to save data
handles.neighbors = trace_curr.neighbors;

% set all the things normally set in loading a new data file
handles.xydisplayed = 1;
if isfield(trace_curr,'time')
handles.time = trace_curr.time;
end
set(handles.FileName,'string',filename);
set(handles.FilePath,'string',path);
set(handles.StepsFilename,'string',strcat(path,filename));