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
step_thresh = 7.5;
avg_window = 12; %8
minsteplength = 5; %1
maxnanratio = 0.5;

% % good params for high temporal resolution dynein MINFLUX
% step_thresh = 7;
% avg_window = 30;
% minsteplength = round(7.0/0.33);
% maxnanratio = 0.5;
% 
% % good params for high temporal resolution 2C dynein MINFLUX 12/03
% step_thresh = 7;
% avg_window = 15;
% minsteplength = round(7.0/0.67);
% maxnanratio = 0.75;

% good params for high temporal resolution 2C kinesin MINFLUX 25/01/01
step_thresh = 4;
avg_window = 12; % accounts for NaNs, helpful if a reasonable size of the step
minsteplength = 12;
maxnanratio = 0.75;

% good params for high temporal resolution 2C kinesin MINFLUX 25/01/01
step_thresh = 5.5;
avg_window = 15; % accounts for NaNs, helpful if a reasonable size of the step
minsteplength = 15;
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
stepidx = find(handles.stepVector(idx(2:end)) - handles.stepVector(idx(1:end-1)));
% idx(stepidx)
stepidx = [1, stepidx, length(idx)];

% % first remove steps that are too short
% for i = length(idx):-1:2
%     % if ~isnan(handles.stepVector(i-1)) && ~isnan(handles.stepVector(i))
%         if abs(handles.stepVector(idx(i)) - handles.stepVector(idx(i-1))) > 0
%             [idx(i-1),idx(i)]
%             % if the number of points in the step is less than a certain
%             % number, remove the one that leads to a smaller change
%             if sum(abs( handles.stepVector(idx(i):min(idx(i-1)+minsteplength, length(handles.stepVector))) - handles.stepVector(idx(i))),'omitnan') > 1e-30
%                 % fprintf("PHASE 1 \n")
%                 % idx(i)
%                 if abs(handles.stepVector(idx(i)) - handles.stepVector(idx(i)-1)) < abs(handles.stepVector(idx(i-1)+1) - handles.stepVector(idx(i-1)))
%                     handles = AddRmvStepManually(hObject, handles, handles.t(idx(i)));
%                 else
%                     handles = AddRmvStepManually(hObject, handles, handles.t(idx(i-1)));
%                 end
%             end
%         end
%     % end
% end


% JS Edit 2025/01/01
% remove the last step if too short
if idx(stepidx(end)) - idx(stepidx(end-1)) < minsteplength
    handles = AddRmvStepManually(hObject, handles, handles.t(idx(stepidx(end-1))));
end
% first remove steps that are too short
for i = length(stepidx)-1:-1:3 %ignore the ends
    % if the number of points in the step is less than a certain
    % number, remove the one that leads to a smaller change
    if idx(stepidx(i)) - idx(stepidx(i-1)) < minsteplength
        % fprintf("PHASE 1 \n")
        % idx(stepidx(i))
        % handles.t(idx(stepidx(i+1)))
        % handles.stepVector(idx(stepidx(i+1)))
        % handles.t(idx(stepidx(i)))
        % handles.stepVector(idx(stepidx(i)))
        % handles.t(idx(stepidx(i-1)))
        % handles.stepVector(idx(stepidx(i-1)))
        if abs(handles.stepVector(idx(stepidx(i+1))) - handles.stepVector(idx(stepidx(i)))) < abs(handles.stepVector(idx(stepidx(i))) - handles.stepVector(idx(stepidx(i-1))))
            % handles.t(idx(stepidx(i)))
            handles = AddRmvStepManually(hObject, handles, handles.t(idx(stepidx(i))));
        else
            % handles.t(idx(stepidx(i)))
            handles = AddRmvStepManually(hObject, handles, handles.t(idx(stepidx(i-1))));
        end
    end
end
% remove the first step if too short
if idx(stepidx(2)) - idx(stepidx(1)) < minsteplength
    handles = AddRmvStepManually(hObject, handles, handles.t(idx(stepidx(2))));
end

% stepidx = find(handles.stepVector(idx(2:end)) - handles.stepVector(idx(1:end-1)));
% idx(stepidx)

% JS Edit 2025/01/01 works for 2C now too
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
            
            % We are going to account for NaNs in the average window!!
            

            % find mean to remove
            % fprintf(strcat(num2str(i),'\n'))
            % 
            % idx(max([i-1-avg_window,1]))
            % idx(i-1)
            % idx(i)
            % idx(min([i+avg_window,length(idx)]))

            % if abs( mean(PSD1Data_LongPadded(idx(i):idx(i-1+avg_window)),'omitnan') - mean(PSD1Data_LongPadded(idx(i+avg_window):idx(i-1+2*avg_window)),'omitnan') ) < step_thresh
            if abs( mean(PSD1Data_LongPadded(idx(max([i-1-avg_window,1])):idx(i-1)),'omitnan') - mean(PSD1Data_LongPadded(idx(i):idx(min([i+avg_window,length(idx)]))),'omitnan') ) < step_thresh
                % Then remove said steps
%                 fprintf("Psych no step \n")
                % fprintf("PHASE 2 \n")
                % idx(i)
                handles = AddRmvStepManually(hObject, handles, handles.t(idx(i)-1));
            end

            % if >75% of data points in the region around the step are NaN
            % values, remove because low data quality
            % if sum(isnan(PSD1Data_LongPadded(idx(i):idx(i)-1+avg_window))) > maxnanratio*avg_window || sum(isnan(sum(isnan(PSD1Data_LongPadded(idx(i):idx(i)-1+avg_window))))) > maxnanratio*avg_window
            if sum(isnan(PSD1Data_LongPadded(idx(max([i-1-avg_window,1])):idx(i-1)))) > maxnanratio*length(idx(max([i-1-avg_window,1])):idx(i-1)) || sum(isnan(PSD1Data_LongPadded(idx(i):idx(min([i+avg_window,length(idx)]))) ) ) > maxnanratio*length(idx(i):idx(min([i+avg_window,length(idx)])))
                % Then remove said steps
%                 fprintf("Psych no step \n")
                % fprintf("PHASE 3 \n")
                % idx(i)
                handles = AddRmvStepManually(hObject, handles, handles.t(idx(i)-1));
            end

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

