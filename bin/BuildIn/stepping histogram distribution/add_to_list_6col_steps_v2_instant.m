function [x_steps, y_steps] = add_to_list_6col_steps_v2_instant(varargin)

% JS adaptation of add_to_list_6col_steps_v2 so that it will instead look
% for an "instantaneous jump magnitude" by averaging only raw data points
% closest to jump. This is especially important for the off-axis.

window_size = 10;

trace = varargin{1};
use_thresh = varargin{2};

if nargin >= 3
    x_steps = varargin{3};
    y_steps = varargin{4};
else
    x_steps = [];
    y_steps = [];
end


limit = length(trace(:,1));
counter = length(x_steps);

for i=2:limit
    if trace(i,5) >= 1 && abs(trace(i,3) - trace(i-1,3)) >= use_thresh
        counter  = counter+1;
        % first check whether the step is longer than window_size before
        % and after. If it isn't then just use the already existing average
        % since that would do the same thing

        if i-window_size < 1 %you hit the beginning of the trace
            beforeavg = trace(i-1,3);
        elseif sum(trace(i-window_size:i-1,5)) > 0
            beforeavg = trace(i-1,3);
        else
            beforeavg = mean(trace(i-window_size:i-1,1),'omitnan');
        end

        if i-1+window_size > limit %you hit the end of the trace
            afteravg = trace(i,3);
        elseif sum(trace(i:i-1+window_size,5)) > 0
            afteravg = trace(i,3);
        else
            afteravg = mean(trace(i:i-1+window_size,1),'omitnan');
        end

        x_steps(counter) = afteravg - beforeavg;
        y_steps(counter) = trace(i,4) - trace(i-1,4); %just ignoring it anyway
    end
end