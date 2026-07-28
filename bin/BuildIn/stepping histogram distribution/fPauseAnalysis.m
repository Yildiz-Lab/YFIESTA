function [normalized_pause_density, pause_fraction, V, pause_events, numside, numback, num_start_trunc, num_end_trunc, distance] = fPauseAnalysis(trace, trace_yx, tchoice, threshold_cnt_pause, region_size)
% function [V, normalized_pause_frequency] = fPauseAnalysis(trace, tchoice, threshold_cnt_pause)

if nargin < 5 || isempty(region_size)
    region_size = Inf; %DC260727 no region context supplied -- always use empirical distance
end

% Similar to add_to_list_6col
fulltrace_length = size(trace,1); %DC260727 full (pre-subset) trace length, for start/end truncation checks
ydata_full = trace(:,1); %DC260727 full (pre-subset) on-axis data, for adjacent out-of-region point lookups
ydata = trace(:,1);
if isempty(tchoice)
    tchoice = 1:length(ydata);
end
[idx,~] = find(~isnan(trace(:,3)));
nntchoice = intersect(tchoice, idx); %get no NaN tchoice
if isempty(nntchoice) %DC260727 guard against a region assignment with zero valid points (e.g. only the trimmed-off last frame)
    normalized_pause_density = NaN;
    pause_fraction = NaN;
    V = [];
    pause_events = 0;
    numside = 0;
    numback = 0;
    num_start_trunc = 0;
    num_end_trunc = 0;
    fprintf("empty")
    distance = NaN;
    return;
end
ydata = ydata(nntchoice);
trace = trace(nntchoice,:);
trace_yx = trace_yx(nntchoice,:);

%Then histogram (in code form) and look for pause events
% binsize = 50; %nm
% binsize = 25;
binsize = 100;
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
%             if sum(trace_yx(idx,5)) > 0 % only one count per any pause
%                 numside = numside + 1;
%             end
%             if any(trace(idx+1,3)-trace(idx,3) < 0 ) % only one count per any pause
%                 numback = numback + 1;
%             end 
            % count each individually
            numside = numside + sum(trace_yx(idx,5)); 
            numback = numback + sum(trace(idx+1,3)-trace(idx,3) < 0);
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
num_start_trunc = 0; %DC260727 count of chunks truncated by trace start
num_end_trunc = 0; %DC260727 count of chunks truncated by trace end
for j=1:length(timejump)-1
    chunkidx = timejump(j):timejump(j+1);
    firstpt = nntchoice(timejump(j));
    lastpt = nntchoice(timejump(j+1));
    is_start_trunc = (firstpt == 1);
    is_end_trunc = (lastpt == fulltrace_length);
    if is_start_trunc
        num_start_trunc = num_start_trunc + 1;
    end
    if is_end_trunc
        num_end_trunc = num_end_trunc + 1;
    end

    if ~is_start_trunc && ~is_end_trunc && isfinite(region_size)
        %DC260727 fully-observed crossing of a bounded region
        adddist = region_size;
    else
        %DC260727 truncated chunk (either or both sides), or unbounded catch-all region:
        % use the observed excursion within the chunk itself first
        adddist = max(ydata(chunkidx),[],'omitnan') - min(ydata(chunkidx),[],'omitnan');
        if isnan(adddist)
            adddist = 0;
        end
        if adddist == 0
            %DC260727 zero-spread chunk: fall back to distance from the connecting
            % out-of-region point, on whichever side isn't trace-truncated
            candidates = [];
            if firstpt > 1
                candidates(end+1) = abs(ydata(chunkidx(1)) - ydata_full(firstpt-1));
            end
            if lastpt < fulltrace_length
                candidates(end+1) = abs(ydata(chunkidx(1)) - ydata_full(lastpt+1));
            end
            if ~isempty(candidates)
                adddist = max(candidates);
            end
            %DC260727 if no candidates exist (both sides truncated AND zero spread),
            % adddist legitimately stays 0 -- a fully stationary trace contained in the region
        end
    end
    distance = distance + adddist;
end
normalized_pause_density = pause_events / distance;

end

