function [normalized_pause_frequency, V] = fPauseAnalysis(trace, tchoice, threshold_cnt_pause)
% function [V, normalized_pause_frequency] = fPauseAnalysis(trace, tchoice, threshold_cnt_pause)

% Similar to add_to_list_6col
ydata = trace(:,1);
if isempty(tchoice)
    tchoice = 1:length(ydata);
end
[idx,~] = find(~isnan(trace(:,3)));
nntchoice = intersect(tchoice, idx); %get no NaN tchoice
ydata = ydata(nntchoice);

%Then histogram (in code form) and look for pause events
binsize = 40; %nm
bins = min(ydata):binsize:max(ydata)+binsize; %+binsize to pick up stragglers in case number is >binsize
V = zeros(1,length(bins)-1); %we fill in between values

for i = 1:length(bins)-1
    yhist = ydata(ydata >= bins(i));
    yhist = yhist(yhist < bins(i+1));
    V(i) = length(yhist);
end

%finally analyze for if the value is above the "pause event" threshold
pause_events = length(V(V > threshold_cnt_pause));
% normalized_pause_frequency = pause_events;
% normalized pause frequency is number of pause events per some number of
% time points. To make the numbers more manageable, we will normalize it to
% number of pause events / __ number of time points
norm_time_pts = 100;
normalized_pause_frequency = pause_events / length(nntchoice) * norm_time_pts;

end

