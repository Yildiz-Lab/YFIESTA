function x3 = StepThresholdRemoval_op(x1,x3,t,st,aw,ml)
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


step_thresh = st;
avg_window = aw;
minsteplength = ml;


if aw > length(x1)/3  % Prevent oversized averaging windows
    aw = floor(length(x1)/3);
    warning('Reduced avg_window to %d for data length %d', aw, length(x1));
end

x1 = horzcat(nan(1,avg_window),x1,nan(1,avg_window));
x3 = horzcat(nan(1,avg_window),x3,nan(1,avg_window));



idx = find(~isnan(x3));
stepidx = find(x3(idx(2:end)) - x3(idx(1:end-1)));
% idx(stepidx)

stepidx = [1, stepidx, length(idx)];



if sum(~isnan(x1)) < 2*aw  % Insufficient non-NaN data
    return;  % Skip processing
end

% JS Edit 2025/01/01
% remove the last step if too short
if idx(stepidx(end)) - idx(stepidx(end-1)) < minsteplength
    handles = AddRmvStepManually_op(x1,x3,idx(stepidx(end-1)));
end
% first remove steps that are too short
for i = length(stepidx)-1:-1:3 %ignore the ends
    
    if idx(stepidx(i)) - idx(stepidx(i-1)) < minsteplength
        if abs(x3(idx(stepidx(i+1))) - x3(idx(stepidx(i)))) < abs(x3(idx(stepidx(i))) - x3(idx(stepidx(i-1))))
            % handles.t(idx(stepidx(i)))
            x3 = AddRmvStepManually_op(x1,x3,idx(stepidx(i)));
        else
            % handles.t(idx(stepidx(i)))
            x3 = AddRmvStepManually_op(x1,x3,idx(stepidx(i-1)));
        end
    end
end
% remove the first step if too short
if idx(stepidx(2)) - idx(stepidx(1)) < minsteplength
    x3 = AddRmvStepManually_op(x1,x3,idx(stepidx(1)));
end

% stepidx = find(x3(idx(2:end)) - x3(idx(1:end-1)));
% idx(stepidx)

% JS Edit 2025/01/01 works for 2C now too
% then remove steps based on size

for step_idx = 2:length(stepidx)-1 
    i_step = stepidx(step_idx); % Position in idx array where step occurs
    
    % Calculate segment indices around the step
    prev_segment_start = max(1, i_step - avg_window);
    prev_segment_end = i_step - 1;
    next_segment_start = i_step;
    next_segment_end = min(length(idx), i_step + avg_window);
    
    % Get original x3 indices (not just idx array positions)
    x3_prev_start = idx(prev_segment_start);
    x3_prev_end = idx(prev_segment_end);
    x3_next_start = idx(next_segment_start);
    x3_next_end = idx(next_segment_end);
    
    % Calculate step size using original data (x1)
    prev_mean = mean(x1(x3_prev_start:x3_prev_end), 'omitnan');
    next_mean = mean(x1(x3_next_start:x3_next_end), 'omitnan');
    step_size = abs(next_mean - prev_mean);
    
    % Threshold check
    if step_size < step_thresh
        % Remove step at the original x3 index
        x3 = AddRmvStepManually_op(x1, x3, idx(i_step));
    end
end
end

