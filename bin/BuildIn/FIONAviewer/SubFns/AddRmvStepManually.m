function app = AddRmvStepManually(app, nearestX, nearestY)
% Adds or removes a step where the user clicked on the axes

app = getappdata(app.fig, 'app');
% JS Edit 2024/12/04 to ignore Nans in stepVector
% Before we were just interpolating. Now we are actually doing work
if strcmp(app.FitPanelFields.DropFitAxis.Value,'Long-axis')
    stepVector = app.Data.stepVector;
    rawData = app.Data.PSD1Data_Long;
else
    stepVector = app.Data.shortstepVector;
    rawData = app.Data.PSD1Data_Short;
end
idx = find(~isnan(stepVector));


%"clicked" index is at indexT
%so either add a step at indexT
%or remove one at the step nearest indexT...right?
[~,indexT] = min(abs(app.Data.t - nearestX));
currY = stepVector(indexT);

% MD: fixed offset problem...and uncoupled this process from fit-marking
% I think the indexing bug is fixed right here with those + and -1 bits.
% find index of previous step:
prevStepIndex = find(stepVector(1:indexT)-currY, 1, 'last') + 1;

% find index of next step:
nextStepIndex = find(stepVector(indexT:end)-currY, 1, 'first') + indexT - 1;

% determine if a step is being added or deleted:
if app.FitPanelFields.BtnAdd.UserData.pressed

    % JS Edit 2025/01/01 so update doesn't change the number of Nan's which
    % throws everything off, just use idx instead

    idxlimit = length(idx);
    % changespots finds spots where idx(2:idxlimit)
    changespots = find(stepVector(idx(2:idxlimit)) - stepVector(idx(1:idxlimit-1)));
    [~, k] = min(abs(idx - indexT));
    changespots = [changespots k];
    changespots = sort(changespots);

    changespots = [0, changespots, idxlimit]; %add beginning and end for padding
    for i = 2:length(changespots) % go through and fill in with the average
        % idx(changespots(i-1)+1:changespots(i))
        % [idx(changespots(i-1)+1), idx(changespots(i))]
        % mean(rawData(idx(changespots(i-1)+1:changespots(i))),'omitnan');
        stepVector(idx(changespots(i-1)+1:changespots(i))) = mean(rawData(idx(changespots(i-1)+1:changespots(i))),'omitnan');

    end
    app = UserStepChangesTable(app, k+1, "user");

    
else % remove a step
    
%       MD: Make a changepoints array, = 1 if theres a step, 0 otherwise. No
%       NaNs can be used to make a changepoint
    if (nextStepIndex-indexT) > (indexT-prevStepIndex)
        remove_index = prevStepIndex;
    else
        remove_index = nextStepIndex;
    end

    % JS Edit 2025/01/01 so update doesn't change the number of Nan's which
    % throws everything off, just use idx instead

    idxlimit = length(idx);
    % changespots finds spots where idx(2:idxlimit)
    changespots = find(stepVector(idx(2:idxlimit)) - stepVector(idx(1:idxlimit-1)));
    [~, k] = min(abs(idx(changespots+1) - remove_index));
    app = UserStepChangesTable(app, changespots(k)+1, "user");
    changespots(k) = [];

    changespots = [0, changespots, idxlimit]; %add beginning and end for padding
    for i = 2:length(changespots) % go through and fill in with the average
        % idx(changespots(i-1)+1:changespots(i))
        % [idx(changespots(i-1)+1), idx(changespots(i))]
        % mean(rawData(idx(changespots(i-1)+1:changespots(i))),'omitnan');
        stepVector(idx(changespots(i-1)+1:changespots(i))) = mean(rawData(idx(changespots(i-1)+1:changespots(i))),'omitnan');

    end
    
end

% Finally update the trace data in package trace
stepFit.StepFit = stepVector;
app = PackageTrace(app, stepFit);
% and package into plotted function
if strcmp(app.FitPanelFields.DropFitAxis.Value,'Long-axis')
    app.Data.stepVector = stepVector;
else
    app.Data.shortstepVector = stepVector;
end

setappdata(app.fig, 'app', app);
