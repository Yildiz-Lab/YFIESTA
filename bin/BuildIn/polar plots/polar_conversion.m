function [r, theta] = polar_conversion(data, plot_individual)

% JS 2024/03/27

%% Description
%     convert on and off axis stepping into an r, theta coordinate system

% Parameters
%     data : 
%     Take in data structure which has subfields, most importantly
%     trace and trace_yx, Nx6 structures where
%     (:,1) major axis data
%     (:,2) minor axis data
%     (:,3) step fit to major axis data
%     (:,4) step fit to minor axis data
%     (:,5) places where major axis changes
%     trace means on-axis is major axis, trace_yx means off-axis is major axis
% 
%     plot_individual : (default false) boolean to plot for every run
% 
%     window : integer averaging window for the case of steps

% Returns
%     [r, theta] : coordinates of r,theta from x,y

%% Code
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

% if want to merge_step_components
data = merge_step_components(data);

idx = find(data.trace_2d(:,5)>0); %idx = find(data.trace(:,5)>0);
r = nan(1,length(idx));
theta = nan(1,length(idx));
idx = [1; idx; length(x)]; %pad so that we can refer to steps i-1 and i+1 near beginning and end respectively

% then we just need to go through each step and ask ourselves what is the
% mean distribution

window = 25;

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