function data =  merge_step_components(data)

% JS 2024/04/05

%% Description
% merge on-off axis steps if they belong to the same event of stepping
% within a certain index window defined by merge_window

% Parameters
% data : 
% Take in data structure which has subfields, most importantly
% trace and trace_yx, Nx6 structures where
% (:,1) major axis data
% (:,2) minor axis data
% (:,3) fit to major axis data
% (:,4) fit to minor axis data
% (:,5) places where major axis changes
% trace means on-axis is major axis, trace_yx means off-axis is major axis

% merge_window : integer for largest separation where can count as same
% step

% Returns
% data : updated data structure with new field trace_2d that has new
% positions and means about common cp

%% Code
% want to create the structure such that steps aren't broken into
% components, an xy total step if you will. We will look at both x and y
% direction.

merge_window = 4; % if step occurs within this span, consider it the same step.

% Find steps from data.trace and data.trace_yx
onaxis_idx = find(data.trace(:,5));
offaxis_idx = find(data.trace_yx(:,5));
% length(offaxis_idx)
for i = 1:length(onaxis_idx)
    [d, idx] = min(abs(onaxis_idx(i) - offaxis_idx));
    if d <= merge_window %erase this index, it isn't uniquely only sideways
        offaxis_idx(idx) = [];
    end
end
% length(offaxis_idx)

% Now go through and set these new steps by calculating average around
% them. Also set a threshold that if the size is too small we should ignore

data.trace_2d = data.trace;
offadd = zeros(length(data.trace),1);
offadd(offaxis_idx) = 1;
data.trace_2d(:,5) = data.trace_2d(:,5) + offadd;

% now we just have to go through wherever there is a one, and mean together
% the steps between ones

allsteps_idx = find(data.trace_2d(:,5));
allsteps_idx = [1; allsteps_idx; length(data.trace_2d)];
for j = 2:length(allsteps_idx)-1
    mask = allsteps_idx(j-1):allsteps_idx(j)-1;
    
    % means
    data.trace_2d(mask,3) = mean(data.trace_2d(mask,1),'omitnan');
    data.trace_2d(mask,4) = mean(data.trace_2d(mask,2),'omitnan');
    
end


