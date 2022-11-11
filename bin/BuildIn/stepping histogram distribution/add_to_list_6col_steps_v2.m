function [x_steps, y_steps] = add_to_list_6col_steps_v2(varargin)

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

for (i=2:limit)
    if ( trace(i,5) >= 1 && prod(trace(i-1:i,6))>=use_thresh^2 )
        counter  = counter+1;
        x_steps(counter) = trace(i,3) - trace(i-1,3);
        y_steps(counter) = trace(i,4) - trace(i-1,4);
    end
end