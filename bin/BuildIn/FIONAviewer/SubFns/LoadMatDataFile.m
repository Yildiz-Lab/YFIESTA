function app = LoadMatDataFile(app, path, filename)
% LoadMatDataFile  (modernized for UIFIGURE app)
%
% app      – struct containing UI handles + application data
% path     – folder path
% filename – file name to load

% Reset step fitter
app.Data.stepVector = [];
app.Data.shortstepVector = [];

%% -----------------------------
% Load MAT file
% ------------------------------
DataFilePath = fullfile(path, filename);
trace_curr = load(DataFilePath);

% New MAT format wrapper
if isfield(trace_curr,'data')
    trace_curr = trace_curr.data;
end
app.data = trace_curr; % just load everything for a template

%% -----------------------------
% Load raw coordinate traces
% ------------------------------
app.Data.PSD1Data_Long  = trace_curr.xy(:,1)';
app.Data.PSD1Data_Short = trace_curr.xy(:,2)';
app.Data.t = 1:length(trace_curr.xy(:,1));
if isfield(trace_curr, 'time')
    app.Data.t = trace_curr.time - trace_curr.time(1);
    app.data.time = trace_curr.time;
    app.time_bool = 1;
end

%% Step/vector data (optional depending on file contents)
if isfield(trace_curr,'trace')
    app.data.trace = trace_curr.trace;
    app.Data.stepVector       = app.data.trace(:,3)';
    app.Data.shortstepVector  = app.data.trace(:,4)';
else
    app.Data.stepVector = [];
    app.Data.shortstepVector = [];
end

%% Neighbors field
if isfield(trace_curr,'neighbors')
    app.data.neighbors = trace_curr.neighbors;
else
    app.data.neighbors = [];
end

%% -----------------------------
% Update UI labels
% ------------------------------
app.FileName.Text     = filename;
app.FilePath.Text     = path;
app.LeftPanelFields.StepsFilename.Text = fullfile(path,filename);

