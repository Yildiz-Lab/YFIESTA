function trace = modify_filter_trace(varargin)

trace = varargin{1};
use_thresh = varargin{2}; %in s
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

% figure()
% hold on
% plot(time,trace(:,1),'k-')
% plot(time,trace(:,3),'b-')

idx = find(trace(:,5)>0);
idx = [1; idx; length(idx)];
for j = 3:length(idx)-1
    if time(idx(j)) - time(idx(j-1)) < use_thresh %we are going to remove this step
        % if abs(trace(idx(j),3) - trace(idx(j-1),3)) < abs(trace(idx(j-1),3) - trace(idx(j-2),3))
            %if previous step is further, then merge with later
            % trace(idx(j-1):idx(j+1)-1,3) = mean(trace(idx(j-1):idx(j+1)-1,1),'omitnan')*ones(idx(j+1)-idx(j-1),1);
            % trace(idx(j-1):idx(j+1)-1,4) = mean(trace(idx(j-1):idx(j+1)-1,2),'omitnan')*ones(idx(j+1)-idx(j-1),1);
            % trace(idx(j),5) = 0;
        % else % for all
        if trace(idx(j-1),3) < trace(idx(j-2),3) && trace(idx(j),3) > trace(idx(j-1),3) %only backsteps that have a forward step after (dip)
            % merge with previous
            % trace(idx(j-2):idx(j)-1,3) = mean(trace(idx(j-2):idx(j)-1,1),'omitnan')*ones(idx(j)-idx(j-2),1);
            % trace(idx(j-2):idx(j)-1,4) = mean(trace(idx(j-2):idx(j)-1,2),'omitnan')*ones(idx(j)-idx(j-2),1);
            % trace(idx(j-1),5) = 0;
            % merge with later
            trace(idx(j-1):idx(j+1)-1,3) = mean(trace(idx(j-1):idx(j+1)-1,1),'omitnan')*ones(idx(j+1)-idx(j-1),1);
            trace(idx(j-1):idx(j+1)-1,4) = mean(trace(idx(j-1):idx(j+1)-1,2),'omitnan')*ones(idx(j+1)-idx(j-1),1);
            trace(idx(j),5) = 0;
        end
    end
end

% plot(time,trace(:,3),'g-')

% we have now adapted trace! Yay!! We can return that.

% start = 1;
% counter = length(dwells);
% for i = 2:limit
%     if (cp(i) == 1)
%         prev_start = start;
%         dwell = prev_start:i-1;
%         start = i;
%         if ( prod(usage(start:i)) >= use_thresh && start ~= 1)
%             counter = counter+1;
%             dwells(counter,1) = trace(i-1,4);
%             dwells(counter,2) = trace(i-1,3);
%             % dwells(counter,3) = length(dwell)*time;
%             dwells(counter,3) = time(i-1)-time(prev_start); %JS Edit 2024/03/07
% 
%             if (trace(i-1,3) < -pi/2 || trace(i-1,3) > pi/2 )
%                 dwells(counter,4) = -1;
%             else
%                 dwells(counter,4) = 1;
%             end
%         end
%     end
% end

end