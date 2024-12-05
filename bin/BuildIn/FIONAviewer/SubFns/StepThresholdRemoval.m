function handles = StepThresholdRemoval(hObject, handles)
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
% step_thresh = 7.5;
% avg_window = 12; %8
% minsteplength = 5; %1
% maxnanratio = 0.5;

% good params for high temporal resolution dynein MINFLUX
% step_thresh = 7;
% avg_window = 30;
% minsteplength = round(7.0/0.33);
% maxnanratio = 0.5;

% good params for high temporal resolution 2C dynein MINFLUX 12/03
step_thresh = 7;
avg_window = 30;
minsteplength = round(7.0/0.33);
maxnanratio = 0.75;

% we have to set this to delete to use the prebuilt AddRmvStepManually func
% remember the old state so we can return it when we plot later
oldstate = get(handles.AddDeleteStep,'String');
set(handles.AddDeleteStep,'String','Delete');

PSD1Data_LongPadded = horzcat(nan(1,avg_window),handles.PSD1Data_Long,nan(1,avg_window));

% JS Edit 2024/12/04
% Remove Nan's entirely from the equation, as there should be no steps
% there anyway.
idx = find(~isnan(handles.stepVector));

% first remove steps that are too short
for i = length(idx):-1:2
    % if ~isnan(handles.stepVector(i-1)) && ~isnan(handles.stepVector(i))
        if abs(handles.stepVector(idx(i)) - handles.stepVector(idx(i-1))) > 0
            % if the number of points in the step is less than a certain
            % number, remove
            if sum(abs( handles.stepVector(idx(i):min(idx(i-1)+minsteplength, length(handles.stepVector))) - handles.stepVector(idx(i))),'omitnan') > 1e-30
                % fprintf("PHASE 1 \n")
                % idx(i)
                handles = AddRmvStepManually(hObject, handles, handles.t(idx(i)));
            end
        end
    % end
end

% then remove steps based on size
for i = 2:length(idx)
    % if ~isnan(handles.stepVector(i-1)) && ~isnan(handles.stepVector(i))
        % check if there is a step
        if abs(handles.stepVector(idx(i)) - handles.stepVector(idx(i-1))) > 0
            
            % check if the average of the points around that step are less
            % than threshold (say window of 5ish unless you are too close to the edge)
            % easier just to pad with NaN for sake of this average
%             handles.t(i-1)
%             abs( mean(PSD1Data_LongPadded(i:i-1+avg_window),'omitnan') - mean(PSD1Data_LongPadded(i+avg_window:i-1+2*avg_window),'omitnan') )
            
            % find mean to remove
            if abs( mean(PSD1Data_LongPadded(idx(i):idx(i)-1+avg_window),'omitnan') - mean(PSD1Data_LongPadded(idx(i)+avg_window:idx(i)-1+2*avg_window),'omitnan') ) < step_thresh
                % Then remove said steps
%                 fprintf("Psych no step \n")
                % fprintf("PHASE 2 \n")
                % idx(i)
                handles = AddRmvStepManually(hObject, handles, handles.t(idx(i-1)));
            end

            % if >75% of data points in the region around the step are NaN
            % values, remove because low data quality
            if sum(isnan(PSD1Data_LongPadded(idx(i):idx(i)-1+avg_window))) > maxnanratio*avg_window || sum(isnan(sum(isnan(PSD1Data_LongPadded(idx(i):idx(i)-1+avg_window))))) > maxnanratio*avg_window
                % Then remove said steps
%                 fprintf("Psych no step \n")
                % fprintf("PHASE 3 \n")
                % idx(i)
                handles = AddRmvStepManually(hObject, handles, handles.t(idx(i-1)));
            end

%             % if >50% of data points in the region around the step are NaN
%             % values, remove because low data quality
%             if sum(isnan(PSD1Data_LongPadded(idx(i):idx(i)-1+avg_window))) > 0.5*avg_window || sum(isnan(sum(isnan(PSD1Data_LongPadded(idx(i):idx(i)-1+avg_window))))) > 0.5*avg_window
%                 % Then remove said steps
% %                 fprintf("Psych no step \n")
%                 handles = AddRmvStepManually(hObject, handles, handles.t(idx(i-1)));
%             end
            
        end
        
    % end
end
% End of JS Edit 2024/12/04

% % first remove steps that are too short
% for i = length(handles.stepVector):-1:2
%     if ~isnan(handles.stepVector(i-1)) && ~isnan(handles.stepVector(i))
%         if abs(handles.stepVector(i) - handles.stepVector(i-1)) > 0
%             % if the number of points in the step is less than a certain
%             % number, remove
%             if sum(abs(handles.stepVector(i:min(i-1+minsteplength,length(handles.stepVector)))-handles.stepVector(i))) > 1e-30
%                 handles = AddRmvStepManually(hObject, handles, handles.t(i));
%             end
%         end
%     end
% end
% 
% % then remove steps based on size
% for i = 2:length(handles.stepVector)
%     if ~isnan(handles.stepVector(i-1)) && ~isnan(handles.stepVector(i))
% 
%         % check if there is a step
%         if abs(handles.stepVector(i) - handles.stepVector(i-1)) > 0
% 
%             % check if the average of the points around that step are less
%             % than threshold (say window of 5ish unless you are too close to the edge)
%             % easier just to pad with NaN for sake of this average
% %             handles.t(i-1)
% %             abs( mean(PSD1Data_LongPadded(i:i-1+avg_window),'omitnan') - mean(PSD1Data_LongPadded(i+avg_window:i-1+2*avg_window),'omitnan') )
% 
%             % find mean to remove
%             if abs( mean(PSD1Data_LongPadded(i:i-1+avg_window),'omitnan') - mean(PSD1Data_LongPadded(i+avg_window:i-1+2*avg_window),'omitnan') ) < step_thresh
%                 % Then remove said steps
% %                 fprintf("Psych no step \n")
%                 handles = AddRmvStepManually(hObject, handles, handles.t(i-1));
%             end
% 
%             % if >50% of data points in the region around the step are NaN
%             % values, remove because low data quality
%             if sum(isnan(PSD1Data_LongPadded(i:i-1+avg_window))) > 0.5*avg_window || sum(isnan(sum(isnan(PSD1Data_LongPadded(i:i-1+avg_window))))) > 0.5*avg_window
%                 % Then remove said steps
% %                 fprintf("Psych no step \n")
%                 handles = AddRmvStepManually(hObject, handles, handles.t(i-1));
%             end
% 
%         end
% 
%     end
% end

set(handles.AddDeleteStep,'String',oldstate);

end

