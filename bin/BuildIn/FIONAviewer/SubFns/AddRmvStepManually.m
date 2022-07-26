function handles = AddRmvStepManually(hObject, handles, clickX)
% Adds or removes a step where the user clicked on the axes

%MD: used my method for keeping track of changepoints
%previous way was hard to follow for me, and was causing weird offsets to
%show up in previous and next steps.
% FilteredPlotOffset = str2double(get(handles.FiltOffsetDistance,'String'));
% if handles.display.filtered
%     dispLongFilt = get(handles.ShowXFilt,'value'); %plot long filtered axis
%     dispShortFilt = get(handles.ShowYFilt,'value'); %plot short filtered axis
% else
%     dispLongFilt = 0;
%     dispShortFilt = 0;
% end
% 
% hSteps = handles.stepLineHandle; % handle to the plotted step line
% hStepsShort = handles.offAxisStepHandle; % plotted short-axis step line
% if handles.FilteredFlag == 1
%     if dispLongFilt == 1
%         hStepsFilt = handles.stepLineHandleFilt;
%     if dispShortFilt == 1
%         hStepsShortFilt = handles.offAxisStepHandleFilt;
% end

scaledStepsT = [1:length(handles.stepVector)]; % in datapoints
%scaledStepsX = get(hSteps, 'YData'); % in nanometers

% find the index of the time point closest to the point the user clicked on
[~, indexT] = min(abs(scaledStepsT - clickX));

%find x-value at current time:
currX = handles.stepVector(indexT);

%"clicked" index is at indexT
%so either add a step at indexT
%or remove one at the step nearest indexT...right?

% MD: fixed offset problem...and uncoupled this process from fit-marking
% I think the indexing bug is fixed right here with those + and -1 bits.
% find index of previous step:
diffPrev = handles.stepVector(1:indexT) - currX;
prevStepIndex = find(diffPrev, 1, 'last') + 1;

% find index of next step:
diffNext = handles.stepVector(indexT:end) - currX;
nextStepIndex = find(diffNext, 1, 'first') + indexT - 1;


% determine if a step is being added or deleted:
addStep = strcmp(get(handles.AddDeleteStep, 'String'), 'Add');
if addStep
    %disp('Adding')
    %this adds a step at indexT...
    
    
    
    %{
    ExpSqNoise = str2double(get(handles.ExpSqNoise, 'String'));
    MinDeltaQ = str2double(get(handles.MinDeltaQ, 'String'));
    MinPtsInStep = str2double(get(handles.MinPtsInStep, 'String'));
        
    % find array to fit steps to:
    x = real(handles.currentPlotPSD_Long(prevStepIndex:nextStepIndex));
    
    % remove NaNs and extrapolate:
    x = RemoveNaNs(x);
    
    MaxNumSteps = 2; % Only fit one step
    % fit the steps
    StepStatistics = ShihBilyardStepFinder_TDV(x,MaxNumSteps,...
        ExpSqNoise,MinDeltaQ,MinPtsInStep,0);
    handles.stepVector(prevStepIndex:nextStepIndex) = ...
        StepStatistics.StepFit;
     %}
   
%     prevStepLevel = ...
%         mean(handles.currentPlotPSD_Long(prevStepIndex:indexT));
%     nextStepLevel = ...
%         mean(handles.currentPlotPSD_Long(indexT:nextStepIndex));


%  MD: Make a changepoints array, = 1 if theres a step, =0 otherwise. No
%  NaNs can be used to make a changepoint
    limit = length(handles.stepVector);
    changepoints = ( (handles.stepVector(2:limit) - handles.stepVector(1:limit-1)) ~= 0 ...
        & isnan(handles.stepVector(2:limit) - handles.stepVector(1:limit-1)) == 0 );
    changepoints = [0 changepoints];
    changepoints(indexT) = 1;
% MD: use the changepoints array to remake the entire on-axis fit...looks
% like that is used right off to make the off-axis fit, so we should solve
% the indexing problem right here. 
% If not, add changepoints to handles, and
% put MD's fit maker in there as well...
    last = 1;
    for (i = 1:limit)
         if (changepoints(i) == 1)
             handles.stepVector(last:i-1) = nanmean(handles.currentPlotPSD_Long(last:i-1));
             last = i;
         end
         if (i == limit)
             handles.stepVector(last:limit) = nanmean(handles.currentPlotPSD_Long(last:limit));
         end
    end
    
    
    % adjust the step array
%     handles.stepVector(prevStepIndex:indexT) = ...
%         ones(1,indexT - prevStepIndex + 1) * prevStepLevel;
%     handles.stepVector(indexT:nextStepIndex) = ...
%         ones(1,nextStepIndex - indexT + 1) * nextStepLevel;
%    
    
else % remove a step
    %disp('Removing')
    
%     if (nextStepIndex-indexT) > (indexT-prevStepIndex)
%     % remove previous step
%         
%         % find index of step before the step being deleted:
%         prevStepX = handles.stepVector(prevStepIndex);
%         diffPrev = handles.stepVector(1:prevStepIndex) - prevStepX;
%         prePrevStepIndex = find(diffPrev, 1, 'last');
%     
%         newStepLevel = ...
%         mean(handles.currentPlotPSD_Long(prePrevStepIndex:nextStepIndex));
%         
%         % adjust the step array
%         handles.stepVector(prePrevStepIndex:nextStepIndex) = ...
%             ones(1,nextStepIndex - prePrevStepIndex + 1) * newStepLevel;
%     
%     else
%     % remove next step
%          % find index of step after the step being deleted:
%         nextStepX = handles.stepVector(nextStepIndex);
%         diffNext = handles.stepVector(nextStepIndex:end) - nextStepX;
%         nextNextStepIndex = find(diffNext, 1, 'first') + nextStepIndex;
%     
%         newStepLevel = ...
%         mean(handles.currentPlotPSD_Long(prevStepIndex:nextNextStepIndex));
%         
%         % adjust the step array
%         handles.stepVector(prevStepIndex:nextNextStepIndex) = ...
%             ones(1,nextNextStepIndex - prevStepIndex + 1) * newStepLevel;
%     
%     end
%       MD: Make a changepoints array, = 1 if theres a step, 0 otherwise. No
%       NaNs can be used to make a changepoint
    if (nextStepIndex-indexT) > (indexT-prevStepIndex)
        remove_index = prevStepIndex;
    else
        remove_index = nextStepIndex;
    end
    %disp (remove_index);
    limit = length(handles.stepVector);
    changepoints = ( (handles.stepVector(2:limit) - handles.stepVector(1:limit-1)) ~= 0 ...
        & isnan(handles.stepVector(2:limit) - handles.stepVector(1:limit-1)) == 0 );
    changepoints = [0 changepoints];
    %disp(changepoints(remove_index-1:remove_index+1));
    changepoints(remove_index) = 0;
% MD: use the changepoints array to remake the entire on-axis fit...looks
% like that is used right off to make the on-axis fit, so we should solve
% the indexing problem right here.
    last = 1;
    for (i = 1:limit)
         if (changepoints(i) == 1)
             handles.stepVector(last:i-1) = nanmean(handles.currentPlotPSD_Long(last:i-1));
             last = i;
         end
         if (i == limit)
             handles.stepVector(last:limit) = nanmean(handles.currentPlotPSD_Long(last:limit));
         end
    end
end




% delete(hSteps)
% delete(hStepsShort)
% delete(hStepsShortFilt)
% delete(hStepsFilt)
%will assume this bit works fine...might want to add a handle for
%changepoints though...
handles.shortStepVector = GenerateShortAxisSteps(handles);

handles = PlotData(hObject, handles, 2);

% hold on
% % Re Plot fitted steps:
% 
% 
%     handles.stepLineHandle = plot( ...
%         handles.stepVector, '-', 'LineWidth', 2, 'Color', [0.2 0.9 0]);
%     
%     handles.offAxisStepHandle = plot( ...
%         handles.shortStepVector, '-', 'LineWidth', 2, 'Color', [0.3 0.3 0.9]);
%     if handles.FilteredFlag == 1 && FilteredPlotOffset ~= 0
%         if dispLongFilt == 1
%             handles.stepLineHandle = plot( ...
%                 handles.stepVector+FilteredPlotOffset, '-', 'LineWidth', 2, 'Color', [0.3 0.3 0.9]);
%         end
%         if dispShortFilt == 1
%             handles.offAxisStepHandle = plot( ...
%                 handles.shortStepVector+FilteredPlotOffset, '-', 'LineWidth', 2, 'Color', [0.3 0.3 0.9]);
%     
%         end
%     end
% 
% hold off

