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
    
    % and deltax by mean of blip and the nearby window or next step
%     x(bidx+n:n:xidx-ceil(0.5*n))
%     time(bidx+n:n:xidx-ceil(0.5*n))
    meanblip = mean(x(bidx+n:n:xidx-ceil(0.5*n)),'omitnan');
    
    next_idx = bidx - n*window;
    
    if abs(bidx - xsteps_idx(idx-n)) < window
        next_idx = xsteps_idx(idx-n) + floor(0.5*n);
        % if ahead, subtract 1 since last spot is at n-1.
        % otherwise, keep at 0
    end
    
    % then we do it
%     bidx:-n:next_idx
    meanstep = mean(x(bidx:-n:next_idx));
    
    dx(i) = meanblip - meanstep;
    
end

end
end


end

