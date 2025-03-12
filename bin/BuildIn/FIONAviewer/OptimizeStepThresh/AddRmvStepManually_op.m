function x3 = AddRmvStepManually_op(x1, x3, targetIndex)
% Adds or removes a step at the target index in the step vector x3
% Inputs:
%   x1 - Original data array
%   x3 - Step vector to modify
%   targetIndex - Index in x3 where the step should be removed
% Output:
%   x3 - Updated step vector with the target step removed

% Find non-NaN indices in x3
idx = find(~isnan(x3));
if isempty(idx)
    return; % No steps to process
end

% Determine the closest non-NaN index to the targetIndex
[~, k] = min(abs(idx - targetIndex));
target_in_idx = idx(k);

% Find all step change positions in the non-NaN indices
step_changes = find(diff(x3(idx)) ~= 0);
if isempty(step_changes)
    return; % No steps to remove
end

% Include boundaries to handle edge cases
step_changes = [0; step_changes(:); length(idx)];

% Find the step change closest to the target_in_idx within the non-NaN indices
[~, closest_step] = min(abs(idx(step_changes(2:end-1) + 1) - target_in_idx));
closest_step = closest_step + 1; % Adjust for the added boundary



% Determine the segments to merge
start_segment = step_changes(closest_step - 1) + 1;
end_segment = step_changes(closest_step + 1);


start_merge = idx(start_segment);
end_merge = idx(end_segment);



if all(isnan(x1(start_merge:end_merge)))
    return;  % Don't modify if segment is all NaN
end


valid_segment = x1(start_merge:end_merge);
valid_segment = valid_segment(~isnan(valid_segment));

if isempty(valid_segment)
    return;  % Don't modify x3 if segment is all NaNs
end



merged_mean = mean(valid_segment);


if isnan(merged_mean)
    merged_mean = mean(x1(start_merge:end_merge), 'omitnan');  % Second attempt
    if isnan(merged_mean)
        return;  % Abort if still NaN
    end
end
% Update the step vector with the merged mean
x3(start_merge:end_merge) = merged_mean;

end
