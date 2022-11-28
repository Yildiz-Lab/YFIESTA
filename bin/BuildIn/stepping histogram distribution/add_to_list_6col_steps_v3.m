function [x_steps, y_steps] = add_to_list_6col_steps_v3(trace,tchoice,use_thresh)
% adapted from add_to_list_6col_steps_v2
x_steps = [];
y_steps = [];

limit = length(trace(:,1));
counter = length(x_steps);

for i=2:limit
    if trace(i,5) >= 1 && abs(trace(i,3) - trace(i-1,3)) >= use_thresh %check that it is a changepoint
        %check that the beginning point is within that region
        % if we want to do end, then check i rather than i-1
        if ismember(i-1,tchoice) 
            counter  = counter+1;
            x_steps(counter) = trace(i,3) - trace(i-1,3);
            y_steps(counter) = trace(i,4) - trace(i-1,4);
        end
    end
end

