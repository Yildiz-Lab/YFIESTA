function [r, theta] = polar_conversion(data, plot_individual)

% data.trace from FIONAviewer x, y, xstep, ystep, step_bool

if nargin < 2
    plot_individual = 0;
end

% I am trying to see the polar plot where we have the forward direction
% with a magnitude angle plot.
% 
% I will just compile the r, theta, then make it into a histogram later

x = data.trace(:,1);
y = data.trace(:,2);

if isfield(data,'trace_yx')
    yidx = find(data.trace_yx(:,5)>0);
    yidx = [1; yidx; length(y)];
else
    yidx = [];
end

idx = find(data.trace(:,5)>0);
r = nan(1,length(idx));
theta = nan(1,length(idx));
idx = [1; idx; length(x)]; %pad so that we can refer to steps i-1 and i+1 near beginning and end respectively

% then we just need to go through each step and ask ourselves what is the
% mean distribution

window = 10;


for i = 2:length(idx)-1
    % i is where the step is, so now we just have to make sure no other
    % steps happen before, we use this window, we're done
    bidx = idx(i)-window-1;
    eidx = idx(i)+window;

    if bidx < idx(i-1)
        bidx = idx(i-1);
    end

    if eidx > idx(i+1)-1
        eidx = idx(i+1)-1;
    end

    % % option to use y trace steps to set end points also
    % if ~isempty(yidx)
    %     % find yidx that could be in this range
    %     find(yidx-window)
    % end

    % Now we just step through calculating <r>, <theta> by window average
    xi = mean(x(bidx:idx(i)-1),'omitnan');
    xf = mean(x(idx(i):eidx),'omitnan');
    yi = mean(y(bidx:idx(i)-1),'omitnan');
    yf = mean(y(idx(i):eidx),'omitnan');
    r(i) = sqrt((xf - xi)^2 + (yf - yi)^2);
    theta(i) = atan((yf-yi)/(xf-xi)) - pi*floor((sign(xf - xi)-1)/2);

    % and repeat forever


end

if plot_individual
plot_polar_conversion(r, theta)
end

end