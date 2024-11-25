function WriteStepsToFile(handles, fileName)
% Creates a new file or appends data to an already open file. Step data is
% stored in the structure stepResults, which contains the following fields:
%
% timeVector: contains the time axis, in seconds
% dataVector: contains the data that the steps are being fit to
% stepVector: contains the fitted steps. Has the same number of points as
%               dataVector and timeVector.
% rss: residual sum of squares; used to determine goodness of fit.
%       Calculated by taking sum over i for i=1:n of (yi-f(xi))^2, where
%       yi = actual data point and f(xi) = value of fitted step at that
%       time point
% rssNorm: same as above, but normalized by number of samples; in other
%       words, rss per data point. rssNorm = rss/n 
% filename: name of the original data file
% runTime: [tBegin, tEnd] - in seconds, the start and end points of the
%       step-fitter's run. Can be used to screen for overlapping regions
% stepTimes: vector containing the start times for each fitted step, in 
%       seconds
% stepSizes: vector containing the sizes, in nanometers, of all steps.
%       Contains th same number of elements as stepTimes, and stepSizes(i)
%       corresponds to stepTimes(i)
% timeStamp: time of fit as a date vector ([year month day hour minute
%       seconds]. Can be used to determine which fits are newer and which
%       are older.

% Created by Vladislav Belyy
% Last updated on 01/03/2011


% Determine if the file already exists:

%MD rewrite 030812
%in the format of my trace object
%6 column matrix: datax datay fitx fity
%for y: looks like he stores short data in handles.shortAxisData
%and handles.shortStepVector
% loadFail = 0;
% 
% try
%     load(fileName);    
% catch %#ok<CTCH>
%     loadFail = 1;
%     disp('This file does not exist; creating a new one');
% end
% 
% if loadFail % generate new stepResults structure
%     stepResults = struct('timeVector', [], 'dataVector', [], ... 
%         'stepVector', [], 'filename', '', ...
%         'shortStepVector',[],'shortDataVector',[],'changePoints',[],'usage',[]);
%     
%     ind = 1; % start writing to the first field of the srtucture
% 
% else
%     ind = length(stepResults)+1; %#ok<NODEF> % add data to stepResults
% end

% JS Edit 2022/11/08 Comment Out because it almost always goes forward
% % Determine which way is up and which way is down
% button = questdlg( ...
%     'In this trace, which direction is "forwards" for the motor?', ...
%     'Select directionality','Up','Down','Up');
% if strcmp(button, 'Up')
%     flipSteps = 1; % do not flip steps
% else
%     flipSteps = -1; % flip steps
% end
flipSteps = 1; % do not flip steps

% Write data:

% iStart = find(noNaNSteps, 1, 'first'); % index of first sample
% iEnd = find(noNaNSteps, 1, 'last'); % index of last sample
% iStart = 1;
% iEnd = length(handles.stepVector);


% stepResults(ind).timeVector = handles.currentPlotT(iStart:iEnd);
% stepResults(ind).dataVector = handles.currentPlotPSD_Long(iStart:iEnd);
% stepResults(ind).stepVector = handles.stepVector(iStart:iEnd);
% 
% stepResults(ind).shortStepVector = handles.shortStepVector(iStart:iEnd);
% stepResults(ind).shortDataVector = handles.currentPlotPSD_Short(iStart:iEnd);
% stepResults(ind).rss = sum((stepResults(ind).dataVector - ...
%     stepResults(ind).stepVector).^2);
% stepResults(ind).rssNorm = stepResults(ind).rss / ...
%                                 length(stepResults(ind).timeVector);


% stepResults(ind).runTime = [stepResults(ind).timeVector(1), ...
%     stepResults(ind).timeVector(end)];

% Calculate stepTimes and stepSizes:

% stepIndices = logical(stepResults(ind).stepVector(2:end) - ...
%     stepResults(ind).stepVector(1:end-1));
% 
% stepResults(ind).stepTimes = stepResults(ind).timeVector(stepIndices);
% stepResults(ind).stepSizes = (stepResults(ind).stepVector([false,...
%     stepIndices])-stepResults(ind).stepVector(stepIndices))*flipSteps;


% stepResults(ind).timeStamp = clock; %#ok<NASGU>

%make "changepoints" and "usage" columns
    limit = length(handles.stepVector);
    changePoints = ( (handles.stepVector(2:limit) - handles.stepVector(1:limit-1)) ~= 0 ...
        & isnan(handles.stepVector(2:limit) - handles.stepVector(1:limit-1)) == 0 );
    changePoints = [0 changePoints];
    
% save data 
%save(fileName, 'stepResults');
trace = [(handles.currentPlotPSD_Long*flipSteps)' (handles.currentPlotPSD_Short*flipSteps)' ...
    (handles.stepVector*flipSteps)' (handles.shortStepVector*flipSteps)'...
    (changePoints)' (handles.usageVector)' ];

% JS Edit 2022/11/11 for .mat files
if handles.xydisplayed
    data.xy = [(handles.currentPlotPSD_Long*flipSteps)' (handles.currentPlotPSD_Short*flipSteps)'];
    data.yx = data.xy(:,[2,1]);
    data.trace = trace;
else
    data.yx = [(handles.currentPlotPSD_Long*flipSteps)' (handles.currentPlotPSD_Short*flipSteps)'];
    data.xy = data.yx(:,[2,1]);
    data.trace_yx = trace;
end

scavenge = load(fileName);
if isfield(scavenge,'data')
    if isfield(scavenge.data, 'trace') && ~handles.xydisplayed
        data.trace = scavenge.data.trace;
    end
    if isfield(scavenge.data, 'trace_yx') && handles.xydisplayed
        data.trace_yx = scavenge.data.trace_yx;
    end
elseif isfield(scavenge,'trace')
    if isfield(scavenge, 'trace') && ~handles.xydisplayed
        data.trace = scavenge.trace;
    end
    if isfield(scavenge, 'trace_yx') && handles.xydisplayed
        data.trace_yx = scavenge.trace_yx;
    end
end

data.neighbors = handles.neighbors;
data.time = handles.time;
%re-worked to get rid of all that stepResults struct business
%this will now re-write traces so be careful!
disp([ 'writing current trace to ' fileName ])
save(fileName,'data');
    