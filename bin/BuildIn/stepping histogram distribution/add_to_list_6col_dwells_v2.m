function dwells = add_to_list_6col_dwells_v2(varargin)

trace = varargin{1};
use_thresh = varargin{2};
dwells = varargin{3};
time = varargin{4};
% if length(time) < 2
%     time = varargin{4}*0:length(trace(:,3))-1;
% end
if length(time) < length(trace(:,3))
    time = (time(2)-time(1))*0:length(trace(:,3))-1;
end
cp = trace(:,5);
usage = trace(:,6);
limit = length(trace(:,3));
%limit = length(trace(:,7)); % JS Edit 220207

start = 1;
counter = length(dwells);
for i = 2:limit
    if (cp(i) == 1)
        prev_start = start;
        dwell = prev_start:i-1;
        start = i;
        if ( prod(usage(start:i)) >= use_thresh && start ~= 1)
            counter = counter+1;
            dwells(counter,1) = trace(i-1,4);
            dwells(counter,2) = trace(i-1,3);
            % dwells(counter,3) = length(dwell)*time;
            dwells(counter,3) = time(i-1)-time(prev_start); %JS Edit 2024/03/07

            if (trace(i-1,3) < -pi/2 || trace(i-1,3) > pi/2 )
                dwells(counter,4) = -1;
            else
                dwells(counter,4) = 1;
            end
        end
    end
end

end