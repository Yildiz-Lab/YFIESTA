function dwells = add_to_list_6col_dwells_v3(trace,tchoice,framerate,use_thresh)
% adapted from add_to_list_6col_dwells_v2
dwells = [];

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
            %check that the beginning point is within that region
            % if we want to do end, then check i rather than i-1
            if ismember(i-1,tchoice) 
                counter = counter+1;
                dwells(counter,1) = trace(i-1,4);
                dwells(counter,2) = trace(i-1,3);
                dwells(counter,3) = length(dwell)*framerate;

                if (trace(i-1,3) < -pi/2 || trace(i-1,3) > pi/2 )
                    dwells(counter,4) = -1;
                else
                    dwells(counter,4) = 1;
                end
            end
        end
    end
end
