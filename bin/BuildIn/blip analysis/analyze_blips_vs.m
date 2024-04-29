function [blipin, blipout, totpts, totsteps] = analyze_blips_vs(data)

% this is basically just find_blips again, but used to report blips that
% were included in the before steps and ones not for sake of p-value

% data struct from fiona and then analyzed by picking blips

window = 30;

% extract the important information

time = data.time;
x = data.trace(:,1);
xsteps = data.trace(:,3);
xsteps_bool = data.trace(:,5);
totsteps = sum(xsteps_bool);

dt = []; dx = []; step = [];

if isfield(data,'blips')
if ~isempty(data.blips)
blip_idx = data.blips(:,1);

xsteps_idx = find(xsteps_bool);
% add in the beginning and end so that we don't ever reach the end of steps
xsteps_idx = [1; xsteps_idx; length(xsteps)];

dt = nan(1,length(blip_idx));
dx = nan(1,length(blip_idx));
step = nan(1,length(blip_idx));

% find step nearest the blip
for i = 1:length(blip_idx)
    bidx = blip_idx(i);
    [~,idx] = min(abs(xsteps_idx - bidx));
    xidx = xsteps_idx(idx);
    
    step(i) = xsteps(xidx)-xsteps(xidx-1);
    
    n = 1; % the step is after the blip
    if xidx < bidx
        n = -1; % the step is before the blip
    end
    
    % now calculate time
    
    dt(i) = n*(time(xidx-ceil(0.5*n)) - time(bidx+n)); 
    if dt(i) < 0 %ignore if dt is negative, there was an error picking blips
        dt(i) = NaN;
    end
    
    % and deltax by mean of blip and the nearby window or next step
%     x(bidx+n:n:xidx-ceil(0.5*n))
%     time(bidx+n:n:xidx-ceil(0.5*n))
    
    % toggle to only include the dip points due to missed steps. So we
    % already know before blip, now we also check that after blip the mean
    % isn't too large to make it positive (so remove ones that cross over the stepping boundary).
    
    % this is just for Ahmet's point
    xblip = x(bidx+n:n:xidx-ceil(0.5*n));
    xfilt = xblip(xblip < xsteps(bidx));
    % xfilt = x(bidx+n:n:xidx-ceil(0.5*n));
    meanblip = mean(xfilt,'omitnan');
    
    % Option for nearest step in size (i.e. before if jump forward, after
    % if jump backward)
%     next_idx = bidx - n*window;
%     
%     if abs(bidx - xsteps_idx(idx-n)) < window
%         next_idx = xsteps_idx(idx-n) + floor(0.5*n);
%         % if ahead, subtract 1 since last spot is at n-1.
%         % otherwise, keep at 0
%     end
    % option for previous step in time (i.e. before for both forward and
    % backward)
%     next_idx = bidx - window;
%     if abs(bidx - xsteps_idx(idx-1)) < window
%         next_idx = xsteps_idx(idx-1);
%     end
    
    % then we do it
%     bidx:-n:next_idx
    if ~isnan(dt(i)) %only include if dt was not negative
%         meanstep = mean(x(bidx:-n:next_idx));
%         meanstep = xsteps(bidx); %get current step in time
        meanstep = xsteps(xidx-1); % get previous step in time
    
        dx(i) = meanblip - meanstep;
    end
    
end

% blip_idx
% step
% dt
% dx


% calculate std and plot 2std
sigma = std(data.trace(:,1) - data.trace(:,3),'omitnan');

% okay look at all data and find points that are >2*sigma dips
two_sigma_mask = data.trace(:,1) < data.trace(:,3)-2*sigma;
dx_mean_subtract = data.trace(:,1) - data.trace(:,3);
step_idx = find(data.trace(:,5) > 0);
step_sign = sign(data.trace(step_idx,3) - data.trace(step_idx-1,3));

% make a mask that goes window before a forward step and window after a
% backward step
window_mask = zeros(length(data.trace(:,1)),1);

set_window = 3;
for i=1:length(step_idx)
    window = set_window;
    % if forward
    if step_sign(i) > 0
        if i-1 < 1
            window = max(1, step_idx(i)-window) - step_idx(i);
        elseif step_idx(i-1) > step_idx(i)-window
            window = step_idx(i)-1 - step_idx(i-1);
        end
        % size(window_mask(step_idx(i)-window:step_idx(i)-1))
        % size(ones(window,1))
        window_mask(step_idx(i)-window:step_idx(i)-1) = ones(window,1);
    elseif step_sign(i) < 0
        if i+1 > length(step_idx)
            window = min(length(step_idx), step_idx(i)+window) - step_idx(i);
        elseif step_idx(i+1)-1 < step_idx(i)+window-1
            window = step_idx(i+1)-1 - step_idx(i);
        end
        % size(window_mask(step_idx(i):step_idx(i)+window-1))
        % size(ones(window,1))
        window_mask(step_idx(i):step_idx(i)+window-1) = ones(window,1);
    end
end

% now we combine to find where the two_sigma_mask and window_mask give us a
% hit on both ends
intersect_mask = and(two_sigma_mask,window_mask);
outsider_mask = and(two_sigma_mask,~window_mask);

% blipin = sum(intersect_mask);
% blipout = sum(outsider_mask);
totpts = length(data.trace(:,1));

blipin = dx_mean_subtract(intersect_mask);
blipout = dx_mean_subtract(outsider_mask);

end
end


end

