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
numsig1 = 0;
% numsig2 = 0;
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
    
    % dt(i) = n*(time(xidx-ceil(0.5*n)) - time(bidx+n))
    dt(i) = n*(time(xidx) - time(bidx+ceil(0.5*n)));
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
    [h1,~] = ttest(xblip - xsteps(bidx+n:n:xidx-ceil(0.5*n)));
    % [h2,~] = ttest2(xblip,xsteps(bidx+n:n:xidx-ceil(0.5*n)));
    if ~isnan(h1)
    numsig1 = numsig1 + h1;
    end
    % if ~isnan(h2)
    % numsig2 = numsig2 + h2;
    % end
    
    
    next_idx = bidx - n*window;
    
    if abs(bidx - xsteps_idx(idx-n)) < window
        next_idx = xsteps_idx(idx-n) + floor(0.5*n);
        % if ahead, subtract 1 since last spot is at n-1.
        % otherwise, keep at 0
    end
    
    % then we do it
%     bidx:-n:next_idx
    if ~isnan(dt(i)) %only include if dt was not negative
        % meanstep = mean(x(bidx:-n:next_idx));
        meanstep = xsteps(bidx);
    
        dx(i) = meanblip - meanstep;
    end
    
end


fprintf(strcat('(', num2str(length(blip_idx)), ',' , num2str(numsig1), ') / ',num2str(sum(xsteps_bool)),'\n'))
% fprintf(strcat(num2str(numsig2),' / ',num2str(sum(xsteps_bool)),'\n'))


% blip_idx
% step
% dt
% dx

end
end


end

