
function app = ScrollPrevNext(app, source) %#ok

% --- Executes on button press in Load_Previous.
% Originally from FIONAviewer file, small edits
% Loads the previous PSDsignals file in the current directory
[path, file, ~] = fileparts(app.LeftPanelFields.StepsFilename.Text);
files = dir(fullfile(path, '*_fiona.mat'));
fileNames = {files.name};

% Find index of the current file
idx = strcmpi(fileNames, strcat(file,'.mat'));
[~, currentIndex] = find(idx);

% index of previous unless if next is pushed
newIndex = currentIndex - 1 + 2*strcmp(source.Text, 'Next  →');

if newIndex > length(idx)
    disp('You''ve reached the end of the directory');
elseif newIndex < 1
    disp('You''ve reached the beginning of the directory');
else
    file = fileNames{newIndex};
    
    % This is synonymous to onLoad
    % Call your loader
    app = LoadMatDataFile(app, path, file);

    % Update plot on the main axes
    app = updateMainPlot(app);
    setappdata(app.fig, 'app', app);

end
