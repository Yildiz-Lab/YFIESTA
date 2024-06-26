function [blipin, blipout, ptsin, ptsout, totpts, totsteps] = analyze_blips_vs_multiple(data,set_window,num2sigma,interval)

% this is basically just find_blips again, but used to report blips that
% were included in the before steps and ones not for sake of p-value

% If you set_window = interval, this is like calculating the number of
% 'steps' that have a preluded 'dip'

% set_window is how far before a forward step or after a backward step we
% look and call a "step"

% this addition asks for the number of 2sigma (num2sigma) points in a
% certain interval (interval) for correlation reporting
% becomes the same as regular vs if num2sigma = 1 and interval = 1

% data struct from fiona and then analyzed by picking blips

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

% set_window = 20;
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
% but now we need to loop over the interval

% % for debugging, use a subset of D1_20366
% debug = 4554:4579;
% find(two_sigma_mask(debug))
% interval_sigma_criteria = zeros(1,length(debug)-interval);
% interval_window_mask = zeros(1,length(debug)-interval);
% outside_window_mask = zeros(1,length(debug)-interval);
% for j = 1:length(debug)-interval+1
%     interval_mask = zeros(length(debug),1);
%     interval_mask(j:j+interval-1) = ones(1,interval);
%     if sum(and(two_sigma_mask(debug),interval_mask)) > 0
%         j
%         sum(and(two_sigma_mask(debug),interval_mask))
%     end
%     interval_sigma_criteria(j) = sum(and(two_sigma_mask(debug),interval_mask)) > num2sigma-1;
%     interval_window_mask(j) = sum(and(interval_mask,window_mask(debug))) > interval-1;
%     outside_window_mask(j) = sum(and(interval_mask,~window_mask(debug))) > interval-1;
% end
% data.trace(debug,5)'
% intersect_mask = and(interval_sigma_criteria,interval_window_mask)
% % outsider_mask = and(interval_sigma_criteria,~interval_window_mask)
% outsider_mask = and(interval_sigma_criteria,outside_window_mask)

interval_sigma_criteria = zeros(1,length(window_mask)-interval);
interval_window_mask = zeros(1,length(window_mask)-interval);
outside_window_mask = zeros(1,length(window_mask)-interval);
for j = 1:length(window_mask)-interval+1
    interval_mask = zeros(length(two_sigma_mask),1);
    interval_mask(j:j+interval-1) = ones(1,interval);
    interval_sigma_criteria(j) = sum(and(two_sigma_mask,interval_mask)) > num2sigma-1;
    interval_window_mask(j) = sum(and(interval_mask,window_mask)) > interval-1;
    outside_window_mask(j) = sum(and(interval_mask,~window_mask)) > interval-1;
end
intersect_mask = and(interval_sigma_criteria,interval_window_mask);
% outsider_mask = and(interval_sigma_criteria,~interval_window_mask);
outsider_mask = and(interval_sigma_criteria,outside_window_mask);

% blipin = sum(intersect_mask);
% blipout = sum(outsider_mask);
% ptsin = sum(interval_window_mask);
% ptsout = length(data.trace(:,1))-ptsin-interval;

totpts = length(data.trace(:,1));

blipin = dx_mean_subtract(intersect_mask);
blipout = dx_mean_subtract(outsider_mask);

ptsin = dx_mean_subtract(and(interval_window_mask,interval_window_mask));
ptsout = dx_mean_subtract(and(outside_window_mask,outside_window_mask));

end
end


end

