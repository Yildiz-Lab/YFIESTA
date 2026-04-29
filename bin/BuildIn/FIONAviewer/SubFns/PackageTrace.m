function app = PackageTrace(app, stepFit)

fill_short_axis = strcmp(app.FitPanelFields.DropFitAxis.Value,'Short-axis');
if fill_short_axis
    % lineObj = findobj(app.AxMain, 'Type', 'Line', 'Tag', 'hShortAxis');
    % lineObjOther = findobj(app.AxMain, 'Type', 'Line', 'Tag', 'hLongAxis');
    y = app.Data.PSD1Data_Short;
    y_other = app.Data.PSD1Data_Long;
else
    % lineObj = findobj(app.AxMain, 'Type', 'Line', 'Tag', 'hLongAxis');
    % lineObjOther = findobj(app.AxMain, 'Type', 'Line', 'Tag', 'hShortAxis');
    y = app.Data.PSD1Data_Long;
    y_other = app.Data.PSD1Data_Short;
end

changepoints_x = [ 0  ( ( stepFit.StepFit(1:length(stepFit.StepFit)-1) - stepFit.StepFit(2:length(stepFit.StepFit)) ) ~= 0) ]; %make a changepoints array (1=step, 0=no step)
% y = lineObj.YData;
if isfield(app.data, 'trace')
    trace = app.data.trace;
else
    trace = NaN(length(y),6);
end

% trace will be set up so that odd numbers are Long-axis and evens are
% Short-axis
% trace_yx is just the odd numbers, just the short events
% col_to_fill is either [1,2,3,5] or [2,1,4,6]

col_to_fill = [1+fill_short_axis, 2-fill_short_axis, 3+fill_short_axis, 5+fill_short_axis];
trace(1:length(y), col_to_fill) = [y' y_other' stepFit.StepFit' changepoints_x' ];

if fill_short_axis
    app.data.trace_yx = NaN(length(y),6);
    app.data.trace_yx(1:length(y), [1 3 5]) = [y' stepFit.StepFit' changepoints_x'];
end
stepidx = find(changepoints_x > 0);

app.data.trace = trace;
app.data.UserChanges = app.Data.StepChangesTable;

end