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
nidx = find(~isnan(trace_curr.xy(:,1)));
app.Data.PSD1Data_Long  = trace_curr.xy(nidx,1)';
app.Data.PSD1Data_Short = trace_curr.xy(nidx,2)';
app.Data.t = 1:length(trace_curr.xy(nidx,1));
if isfield(trace_curr, 'time')
    app.Data.t = trace_curr.time - trace_curr.time(1);
    app.data.time = trace_curr.time;
    app.time_bool = 1;
end

%% Step/vector data (optional depending on file contents)
if isfield(trace_curr,'trace')
    if length(app.data.trace(:,3)) ~= length(app.data.xy(:,1))
        app.data.trace = trace_curr.trace;
    else
        app.data.trace = trace_curr.trace(nidx,:);
    end
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

%% User Changes Table
if isfield(trace_curr,'UserChanges')
    app.Data.StepChangesTable = trace_curr.UserChanges;
else
    app.Data.StepChangesTable = [];
end

%% -----------------------------
% Update UI labels
% ------------------------------
app.FileName.Text     = filename;
app.FilePath.Text     = path;
app.LeftPanelFields.StepsFilename.Text = fullfile(path,filename);

