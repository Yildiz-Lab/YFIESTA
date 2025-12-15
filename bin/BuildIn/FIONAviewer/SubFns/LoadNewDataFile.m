function LoadNewDataFile(app, path, filename)
% Reads data from a new file and updates the app object
%
% Converted for App Designer/uifigure compatibility (2023+)

    %----------------------------------------
    % Turn off filtering
    %----------------------------------------
    app.FilterDataButton.Text = "Filter/Decimate Data?";
    app.display.filtered = 0;

    % Reset step fitter
    app = ResetStepFitter(app);

    % Reset any line-fitting state if applicable
    % app.lineVector = 0;

    % Save full data path
    DataFilePath = fullfile(path, filename);
    app.currentPath = DataFilePath;


    %----------------------------------------
    % Read the data file (your custom function)
    %----------------------------------------
    FIONAData = ReadMatDataFile(DataFilePath);

    if app.xydisplayed
        Data = FIONAData.yx;   % used to be “flip each time”
    else
        Data = FIONAData.xy;
    end

    % Toggle which projection is displayed
    app.xydisplayed = ~app.xydisplayed;

    % Store neighbors & time if present
    if isfield(FIONAData, "neighbors")
        app.neighbors = FIONAData.neighbors;
    else
        app.neighbors = [];
    end

    if isfield(FIONAData, "time")
        app.time = FIONAData.time;
    else
        app.time = [];
    end


    %----------------------------------------
    % Update GUI display fields
    %----------------------------------------
    app.FileNameLabel.Text = filename;
    app.FilePathLabel.Text = path;


    %----------------------------------------
    % Load coordinate data
    %----------------------------------------
    app.PSD1Data_Long  = Data(:,1)';   % long-axis
    app.PSD1Data_Short = Data(:,2)';   % short-axis

    limit = length(Data(:,1));

    % Initialize vectors
    app.stepVector        = nan(1, limit);
    app.shortStepVector   = nan(1, limit);
    app.usageVector       = zeros(1, limit);
    app.changepoints      = zeros(1, limit);

    % Simple index time vector
    app.t = 1:limit;

    % Done
end
