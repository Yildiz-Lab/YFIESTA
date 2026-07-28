function [dt, dx, step, rel_sep_2C] = analyze_blips(data, data2)

% data struct from fiona and then analyzed by picking blips

% Updated to deal with Nan's on 25/03/18 and also two color

window = 15;

% extract the important information
% and remove NaNs if they exist for 2C (blips also assume no NaN)
time = data.time(~isnan(data.time));
x = data.trace(~isnan(data.trace(:,1)),1);
xsteps = data.trace(~isnan(data.trace(:,1)),3);
xsteps_bool = data.trace(~isnan(data.trace(:,1)),5);
if isfield(data,'trace_2d') %we have a merged option so let's use that
    time = data.time(~isnan(data.time));
    x = data.trace_2d(~isnan(data.trace_2d(:,1)), 1);
    xsteps = data.trace_2d(~isnan(data.trace_2d(:,1)), 3);
    xsteps_bool = data.trace_2d(~isnan(data.trace_2d(:,1)), 5);
    ysteps = data.trace_2d(~isnan(data.trace_2d(:,1)), 4);
end

dt = []; dx = []; step = []; rel_sep_2C = [];

if isfield(data,'blips')
if ~isempty(data.blips)
blip_idx = data.blips(:,1);

xsteps_idx = find(xsteps_bool);
% add in the beginning and end so that we don't ever reach the end of steps
xsteps_idx = [1; xsteps_idx; length(xsteps)];

dt = nan(1,length(blip_idx));
dx = nan(1,length(blip_idx));
step = nan(1,length(blip_idx));
rel_sep_2C = nan(1,length(blip_idx));
if nargin > 1 % we have two channels, so we should account for that
    if ~isempty(data2)
    time2 = data2.time(~isnan(data2.time));
    x2 = data2.trace(~isnan(data2.trace(:,1)),1);
    xsteps2 = data2.trace(~isnan(data2.trace(:,1)),3);
    xsteps2_bool = data2.trace(~isnan(data2.trace(:,1)),5);
    if isfield(data2,'trace_2d') %we have a merged option so let's use that
        time2 = data2.time(~isnan(data2.time));
        x2 = data2.trace_2d(~isnan(data2.trace_2d(:,1)), 1);
        xsteps2 = data2.trace_2d(~isnan(data2.trace_2d(:,1)), 3);
        xsteps2_bool = data2.trace_2d(~isnan(data2.trace_2d(:,1)), 5);
        ysteps2 = data2.trace_2d(~isnan(data2.trace_2d(:,1)), 4);
    end
    end
end

% % Option to plot figures and their dips
% figure()
% hold on
% plot(time, x)
% plot(time, xsteps)
% scatter(time(blip_idx), x(blip_idx), 'k', 'filled')
% 
% plot(time2, x2)
% plot(time2, xsteps2)

% find step nearest the blip
numsig1 = 0;
% numsig2 = 0;
for i = 1:length(blip_idx)
    bidx = blip_idx(i);
    [~,bidx_time_nadj] = min(abs(time - data.time(bidx))); % 25/06/20 realized that data.blips refers to the non-nan removed version
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
    % if ~isnan(dt(i)) %only include if dt was not negative
    if dt(i) > 0
%         meanstep = mean(x(bidx:-n:next_idx));
%         meanstep = xsteps(bidx); %get current step in time
        meanstep = xsteps(xidx-1); % get previous step in time
    
        dx(i) = meanblip - meanstep;

        if nargin > 1 % if two channel, include the separation of the meanstep
            if ~isempty(data2)
            % first check that the blip occured while in the intersection
            % of the two data points
            if time(bidx) > min(time2) && time(bidx) < max(time2)
                % find nearest time point in time2 to the initial head
                % separation
                if n > 0
                    [~,idx2] = min(abs(time2 - time(bidx))); %if before, the blip is the start of the step
                    % time(bidx)
                    % time2(idx2)
                else
                    [~,idx2] = min(abs(time2 - time(xidx))); %if after, the step is the start of the step
                    % time(xidx)
                    % time2(idx2)
                end
                % xsteps2(idx2)
                % meanstep
                rel_sep_2C(i) = xsteps2(idx2) - meanstep;
                % Or we could just do it with the merged option
                if isfield(data,'trace_2d')
                    meanstep = [xsteps(xidx-1), ysteps(xidx-1)];
                    rel_sep_2C(i) = sign(rel_sep_2C(i)) * pdist2([xsteps2(idx2), ysteps2(idx2)], meanstep);
                    rel_sep_2C(i) = pdist2([xsteps2(idx2), ysteps2(idx2)], meanstep);
                end
            end
            end
        end
    end
    
end

% time(blip_idx)
% rel_sep_2C

fprintf(strcat('(', num2str(length(blip_idx)), ',' , num2str(numsig1), ') / ',num2str(sum(xsteps_bool)),'\n'))
% fprintf(strcat(num2str(numsig2),' / ',num2str(sum(xsteps_bool)),'\n'))

% blip_idx
% step
% dt
% dx

end
end


end

