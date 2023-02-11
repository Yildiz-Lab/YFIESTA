function [normalized_pause_density, pause_fraction, V, pause_events, numside, numback] = fPauseAnalysis(trace, trace_yx, tchoice, threshold_cnt_pause)
% function [V, normalized_pause_frequency] = fPauseAnalysis(trace, tchoice, threshold_cnt_pause)

% Similar to add_to_list_6col
ydata = trace(:,1);
if isempty(tchoice)
    tchoice = 1:length(ydata);
end
[idx,~] = find(~isnan(trace(:,3)));
nntchoice = intersect(tchoice, idx); %get no NaN tchoice
ydata = ydata(nntchoice);
trace = trace(nntchoice,:);
trace_yx = trace_yx(nntchoice,:);

%Then histogram (in code form) and look for pause events
binsize = 50; %nm
bins = min(ydata):binsize:max(ydata)+binsize; %+binsize to pick up stragglers in case number is >binsize
V = zeros(1,length(bins)-1); %we fill in between values

pause_events = 0;
pause_time = 0;
numback = 0;
numside = 0;
for i = 1:length(bins)-1
    idx = find((ydata >= bins(i)) & (ydata < bins(i+1)));
    idx = idx(idx < length(ydata)); %remove last point so that we don't cause errors searching for steps
    if idx > 0
        %     yhist = ydata(ydata >= bins(i));
        %     yhist = yhist(yhist < bins(i+1));
        V(i) = length(idx);
        if V(i) > threshold_cnt_pause
            pause_events = pause_events + 1;
            pause_time = pause_time + V(i);
            if sum(trace_yx(idx,5)) > 0
                numside = numside + 1;
            end
            if any(trace(idx+1,3)-trace(idx,3) < 0 )
                numback = numback + 1;
            end
        end
    end
end

%finally analyze for if the value is above the "pause event" threshold
% pause_events = length(V(V > threshold_cnt_pause));
% normalized_pause_frequency = pause_events;
% normalized pause frequency is number of pause events per some number of
% time points. To make the numbers more manageable, we will normalize it to
% number of pause events / __ number of time points

norm_time_pts = 100;
normalized_pause_frequency = pause_events / length(nntchoice) * norm_time_pts;
pause_fraction = pause_time / length(nntchoice);

% actually thinking pause frequency should be per run length rather than
% time. Let's see if that speaks to us.
% find where nntchoice jumps
timejump = find(nntchoice(2:end) - nntchoice(1:end-1) > 1)';
timejump = [1, timejump, length(nntchoice)];
distance = 0;
for j=1:length(timejump)-1
    adddist = max(ydata(timejump(j):timejump(j+1))) - min(ydata(timejump(j):timejump(j+1)));
    if isnan(adddist)
        adddist = 0;
    end
    distance = distance + adddist;
end
normalized_pause_density = pause_events / distance;

end

