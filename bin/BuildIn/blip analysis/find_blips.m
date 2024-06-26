function blip_idx = find_blips(fullfilename)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin < 1
    [filename, pathname] = uigetfile('*.mat');
    fullfilename = fullfile(pathname,filename);
end

[d,n,e] = fileparts(fullfilename);
% if contains(f,'._') %JS Edit to delete extra ._ from an error
% f = f(3:end);
% end
savename = fullfile(d,strcat(n,'_wblip',e));

steptrace = load(fullfilename);

data = steptrace.data;

x_line = data.time;
y_line = data.trace(:,1);


% Plot the line
f = figure('KeyPressFcn', @keyPressCallback);
ax = axes('Parent', f);


hold on;
plot(x_line, y_line, 'b.-', 'tag', 'dominant axis');
plot(x_line, data.trace(:,3), 'g-', 'tag', 'dominant fit');

% calculate std and plot 2std
sigma = std(data.trace(:,1) - data.trace(:,3),'omitnan')
plot(x_line, data.trace(:,3)+2*sigma, 'k--', 'tag', '2sigma')
plot(x_line, data.trace(:,3)-2*sigma, 'k--', 'tag', '2sigma')

% Initialize for storage and selecting points
f.UserData.selectPoints = false;
f.UserData.sigma = sigma;
if isfield(data,'blips')
    f.UserData.blips = data.blips;
    % blips already started, appending
    plot(f.UserData.blips(:,2), f.UserData.blips(:,3), 'go', 'MarkerSize', 10);
else
    f.UserData.nearestPoints = [];
    data.blips = [];
end

% okay look at all data and find points that are >2*sigma dips
two_sigma_mask = data.trace(:,1) < data.trace(:,3)-2*sigma;
step_idx = find(data.trace(:,5) > 0);
x_line(step_idx)
step_sign = sign(data.trace(step_idx,3) - data.trace(step_idx-1,3));

% make a mask that goes window before a forward step and window after a
% backward step
window_mask = zeros(length(data.trace(:,1)),1);

set_window = 5; %5 low time res; %20 7.0 power high res; %25 8.0 power high res
for i=1:length(step_idx)
    window = set_window;
    % if forward
    if step_sign(i) > 0
        if i-1 < 1
            window = -max(1, step_idx(i)-window) + step_idx(i);
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
blip_idx_raw = find(intersect_mask > 0);
x_line(two_sigma_mask)
x_line(window_mask > 0)

blip_idx = [];
% finally, we remove extras in a single window spot, we only want the spot
% furthest away
for i = 1:length(step_idx)
    sidx = find(abs(step_idx(i) - blip_idx_raw) < set_window);
    sidx = blip_idx_raw(sidx);
    if ~isempty(sidx)
        if length(sidx) > 1
            if step_sign(i) > 0
                sidx = min(sidx);
            elseif step_sign(i) < 0
                sidx = max(sidx);
            end
        end
        % now find first point that is > -1 sigma (only check one or two
        % points)
        return_sigma = 1.0;
        step_sign(i)
        if step_sign(i) > 0
            bool1 = data.trace(sidx-10:sidx,1) > data.trace(sidx-10:sidx,3)-return_sigma*sigma;
            idx1 = find(bool1 > 0);
            if length(idx1) > 1
                idx1 = idx1(end);
            end
            idx1 = idx1-1+sidx-10;
        elseif step_sign(i) < 0
            bool1 = data.trace(sidx:sidx+10,1) > data.trace(sidx:sidx+10,3)-return_sigma*sigma;
            idx1 = find(bool1 > 0);
            if length(idx1) > 1
                idx1 = idx1(1);
            end
            idx1 = idx1-1+sidx;
        end

    blip_idx = [blip_idx; idx1];
    end

    
end

data.blips = [blip_idx, x_line(blip_idx), y_line(blip_idx)];
blips_percent = length(data.blips)/sum(data.trace(:,5))

save(savename, 'data')

plot(x_line(blip_idx), y_line(blip_idx), 'ko', 'MarkerSize', 50);
f.UserData.blips = data.blips;
f.UserData.blips_percent = size(data.blips,1)/sum(data.trace(:,5));

savefig(fullfile(d,strcat(n,'_blips_sigma',num2str(round(sigma,1)),'.fig')))
close(f)

end

