function [steps, dwells] = two_color_statistics(ch1data,ch2data)
% Take two color data and get statistics:
%   interhead separation and
%   relative rates of stepping
%   step size on interhead separation
%   


fprintf("Now entering the dichromatic regime \n")

% Process the data
if isfield(ch1data, 'time')
    time1 = ch1data.time;
    time1 = ch1data.time(~isnan(ch1data.time));
else
    time1 = 1:size(ch1data.xy,1);
end

figgs = figure();
hold on

idx1 = find(~isnan(ch1data.trace(:,1)));

if isfield(ch2data, 'time')
    time2 = ch2data.time;
    time2 = ch2data.time(~isnan(ch2data.time));
else
    time2 = 1:size(ch2data.xy,1);
end

idx2 = find(~isnan(ch2data.trace(:,1)));
% time1 is not a typo, want to make sure you are normalizing to the
% same time

% On-axis Inter-head separation


% Off-axis Inter-head separation



% Step size


% Dwell



steps = [];
dwells = [];

fprintf("Now leaving the dichromatic regime \n")

end