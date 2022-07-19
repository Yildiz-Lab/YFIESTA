function shortStepVector = GenerateShortAxisSteps(handles)
% Spits out a short axis step vector with step transitions in the same
% places as the long-axis step vector
%
%   Written by Vladislav Belyy on 3-8-2012


stepVector = handles.stepVector;
%shortAxisData = RemoveNaNs(handles.currentPlotPSD_Short);

shortStepVector = stepVector;
%  MD: Make a changepoints array, = 1 if theres a step, =0 otherwise. No
%  NaNs can be used to make a changepoint
    limit = length(handles.stepVector);
    changepoints = ( (handles.stepVector(2:limit) - handles.stepVector(1:limit-1)) ~= 0 ...
        & isnan(handles.stepVector(2:limit) - handles.stepVector(1:limit-1)) == 0 );
    changepoints = [0 changepoints];
% MD: use the changepoints array to remake the entire on-axis fit...looks
% like that is used right off to make the off-axis fit, so we should solve
% the indexing problem right here. 
% If not, add changepoints to handles, and
% put MD's fit maker in there as well...
    last = 1;
    for (i = 1:limit)
         if (changepoints(i) == 1)
             shortStepVector(last:i-1) = nanmean(handles.currentPlotPSD_Short(last:i-1));
             last = i;
         end
         if (i == limit)
             shortStepVector(last:limit) = nanmean(handles.currentPlotPSD_Short(last:limit));
         end
    end

% stepPointsWithNaNs = find(stepVector(2:end) - stepVector(1:end-1));
% 
% stepPoints = stepPointsWithNaNs(~isnan(stepVector(stepPointsWithNaNs))); 
% 
% firstNotNaN = find(~isnan(stepVector), 1, 'first');
% 
% stepPoints = [firstNotNaN stepPoints];
% 
% for i = 2:length(stepPoints)
%     
%     shortStepVector(stepPoints(i-1):stepPoints(i)) = ...
%         mean(shortAxisData(stepPoints(i-1):stepPoints(i)));
% end
% shortStepVector ShortAxisData(
