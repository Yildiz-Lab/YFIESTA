function app = StepThresholdRemoval(app)

%JS 2023/02/10
% Even with twosides code, there is still a lot of overfitting of the
% off-axis data. So the idea is to click a button and apply a manual pure
% "thresholding" to remove steps found that are below this value. This may
% lead to weird mean results, but should get rid of most bogus steps

% Input is same as all subfunctions of FIONAviewer.
% Output should be a new step vector

% User set parameters that should be incorporated into GUI eventually so
% that people can set their own "manual step" corrections

% step_thresh (remove steps smaller than this in the avg_window)
% avg_window (average together points before and after this to determine
% step size)
% minsteplength (minimum number of points in a step to be called a step)

% good params for QD kinesin
step_thresh = 10; %12 %nm
avg_window = 4; %7
minsteplength = 1;
maxnanratio = 0.5;

% good params for LD655 dynein MINFLUX
step_thresh = 7.5;
avg_window = 12; %8
minsteplength = 5; %1
maxnanratio = 0.5;

% good params for high temporal resolution dynein MINFLUX
step_thresh = 7;
avg_window = 30;
minsteplength = round(7.0/0.33);
maxnanratio = 0.5;

% good params for high temporal resolution dynein MINFLUX (by Leo alg)
% step_thresh = 7;
% avg_window = 24;
% minsteplength = 11;
% maxnanratio = 0.5;

% % good params for high temporal resolution 2C dynein MINFLUX 12/03
step_thresh = 7;
avg_window = 15;
minsteplength = round(7.0/0.67);
maxnanratio = 0.75;

% step_thresh = 5;
% avg_window = 10;
% minsteplength = round(7.0);
% maxnanratio = 0.75;

% % good params for high temporal resolution 2C kinesin MINFLUX 25/01/01
% step_thresh = 4;
% avg_window = 12; % accounts for NaNs, helpful if a reasonable size of the step
% minsteplength = 12;
% maxnanratio = 0.75;
% 
% % good params for high temporal resolution 2C kinesin MINFLUX 25/01/01
% step_thresh = 6.5; %5.5
% avg_window = 15; % accounts for NaNs, helpful if a reasonable size of the step
% minsteplength = 15;
% maxnanratio = 0.75;

step_thresh = app.FitPanelFields.MinDeltaChange.Value
avg_window = app.FilterPanelFields.FilterWindow.Value
minsteplength = 5;
maxnanratio = 0.75;



% we have to set this to delete to use the prebuilt AddRmvStepManually func
% remember the old state so we can return it when we plot later
% oldstate = get(handles.AddDeleteStep,'String');
% set(handles.AddDeleteStep,'String','Add');

if strcmp(app.FitPanelFields.DropFitAxis.Value, 'Long-axis')
    stepVector = app.Data.stepVector;
    rawData = app.Data.PSD1Data_Long;
else
    stepVector = app.Data.shortstepVector;
    rawData = app.Data.PSD1Data_Short;
end

% Oh Joseph, you long padded yourself with NaNs
rawDataPadded = horzcat(nan(1,avg_window),rawData,nan(1,avg_window));

% JS Edit 2024/12/04
% Remove Nan's entirely from the equation, as there should be no steps
% there anyway.
idx = find(~isnan(stepVector));
stepidx = find(stepVector(idx(2:end)) - stepVector(idx(1:end-1)));
% idx(stepidx)
stepidx = [1, stepidx, length(idx)];


% JS Edit 2025/01/01
% remove the last step if too short
if idx(stepidx(end)) - idx(stepidx(end-1)) < minsteplength
    app = AddRmvStepManually(app, app.Data.t(idx(stepidx(end-1))), []);
end
% first remove steps that are too short
for i = length(stepidx)-1:-1:3 %ignore the ends
    % if the number of points in the step is less than a certain
    % number, remove the one that leads to a smaller change
    if idx(stepidx(i)) - idx(stepidx(i-1)) < minsteplength
        if abs(stepVector(idx(stepidx(i+1))) - stepVector(idx(stepidx(i)))) < abs(stepVector(idx(stepidx(i))) - stepVector(idx(stepidx(i-1))))
            app = AddRmvStepManually(app, app.Data.t(idx(stepidx(i))),[]);
        else
            app = AddRmvStepManually(app, app.Data.t(idx(stepidx(i-1))),[]);
        end
    end
end
% remove the first step if too short
if idx(stepidx(2)) - idx(stepidx(1)) < minsteplength
    app = AddRmvStepManually(app, app.Data.t(idx(stepidx(2))));
end

% JS Edit 2025/01/01 works for 2C now too
% then remove steps based on size
for i = 2:length(idx)
    % if ~isnan(handles.stepVector(i-1)) && ~isnan(handles.stepVector(i))
        % check if there is a step
        if abs(stepVector(idx(i)) - stepVector(idx(i-1))) > 0

            if abs( mean(rawData(idx(max([i-1-avg_window,1])):idx(i-1)),'omitnan') - mean(rawData(idx(i):idx(min([i+avg_window,length(idx)]))),'omitnan') ) < step_thresh
                disp('steprem')
                app = AddRmvStepManually(app, app.Data.t(idx(i-1)));
            end

            % if >75% of data points in the region around the step are NaN
            % values, remove because low data quality
            % if sum(isnan(PSD1Data_LongPadded(idx(i):idx(i)-1+avg_window))) > maxnanratio*avg_window || sum(isnan(sum(isnan(PSD1Data_LongPadded(idx(i):idx(i)-1+avg_window))))) > maxnanratio*avg_window
            if sum(isnan(rawDataPadded(idx(max([i-1-avg_window,1])):idx(i-1)))) > maxnanratio*length(idx(max([i-1-avg_window,1])):idx(i-1)) || sum(isnan(rawDataPadded(idx(i):idx(min([i+avg_window,length(idx)]))) ) ) > maxnanratio*length(idx(i):idx(min([i+avg_window,length(idx)])))

                app = AddRmvStepManually(app, app.Data.t(idx(i-1)));
            end

        end

    % end
end
% End of JS Edit 2024/12/04

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

end

