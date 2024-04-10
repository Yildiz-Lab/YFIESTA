function [dt, dx, step] = analyze_blips(data)

% data struct from fiona and then analyzed by picking blips

window = 30;

% extract the important information

time = data.time;
x = data.trace(:,1);
xsteps = data.trace(:,3);
xsteps_bool = data.trace(:,5);

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

end
end


end

