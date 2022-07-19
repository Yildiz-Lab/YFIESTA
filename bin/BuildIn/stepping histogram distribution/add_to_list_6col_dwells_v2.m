function dwells = add_to_list_6col_dwells_v2(varargin)

trace = varargin{1};
use_thresh = varargin{2};
dwells = varargin{3};
time = varargin{4};

cp = trace(:,5);
usage = trace(:,6);
limit = length(trace(:,3));
%limit = length(trace(:,7)); % JS Edit 220207

start = 1;
counter = length(dwells);
for i = 2:limit
    if (cp(i) == 1)
        dwell = start:i-1;
        start = i;
        if ( prod(usage(start:i)) >= use_thresh && start ~= 1)
            counter = counter+1;
            dwells(counter,1) = trace(i-1,4);
            dwells(counter,2) = trace(i-1,3);
            dwells(counter,3) = length(dwell)*time;

            if (trace(i-1,3) < -pi/2 || trace(i-1,3) > pi/2 )
                dwells(counter,4) = -1;
            else
                dwells(counter,4) = 1;
            end
        end
    end
end

end