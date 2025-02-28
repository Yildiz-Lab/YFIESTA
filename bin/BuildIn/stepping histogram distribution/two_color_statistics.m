function [steps, dwells] = two_color_statistics(ch1data,ch2data)
% Take two color data and get statistics:
%   interhead separation and
%   relative rates of stepping
%   step size on interhead separation
%   

decimation_factor = 4;

fprintf("Now entering the dichromatic regime \n")

% Process the data
if isfield(ch1data, 'time')
    time1 = ch1data.time;
    % time1 = ch1data.time(~isnan(ch1data.time));
else
    time1 = 1:size(ch1data.xy,1);
end

idx1 = find(~isnan(ch1data.trace(:,1)));
time1 = time1(idx1);


if isfield(ch2data, 'time')
    time2 = ch2data.time;
    % time2 = ch2data.time(~isnan(ch2data.time));
else
    time2 = 1:size(ch2data.xy,1);
end

idx2 = find(~isnan(ch2data.trace(:,1)));
time2 = time2(idx2);


% Ch1 data %x=on-axis  %y=off-axis
time1 = mean_decimate_array(time1, decimation_factor);
x1 = mean_decimate_array(ch1data.trace(idx1,3), decimation_factor);
y1 = mean_decimate_array(ch1data.trace_yx(idx1,3), decimation_factor);
y1 = -y1; %Need to switch to negative due to wrong cross sign for some reason

% Ch1 data changepoints adapted
x1cp = mean_decimate_array(ch1data.trace(idx1,5), decimation_factor); x1cp = ceil(x1cp);
y1cp = mean_decimate_array(ch1data.trace_yx(idx1,5), decimation_factor); y1cp = ceil(y1cp);

% Ch2 data  %x=on-axis  %y=off-axis
time2 = mean_decimate_array(time2, decimation_factor);
x2 = mean_decimate_array(ch2data.trace(idx2,3), decimation_factor);
y2 = mean_decimate_array(ch2data.trace_yx(idx2,3), decimation_factor);

% Ch2 data changepoints adapted
x2cp = mean_decimate_array(ch1data.trace(idx1,5), decimation_factor); x2cp = ceil(x2cp);
y2cp = mean_decimate_array(ch1data.trace_yx(idx1,5), decimation_factor); y2cp = ceil(y2cp);

% Assume Ch1 is shorter
short = 1;
tshort = time1; tlong = time2;
xshort = x1; yshort = y1;
xlong = x2; ylong = y2;
if length(time1) > length(time2) % and then flip if Ch2 is shorter
    short = 2;
    tshort = time2; tlong = time1;
    xlong = x1; ylong = y1;
    xshort = x2; yshort = y2;
end


%% On and off-axis interhead separation

% closest beginning
xdiff = [];
ydiff = [];
for a = 1:length(tshort)
    [~,closestIndex] = min(abs(tlong-tshort(a)));

    % switch whether short is 1 or 2
    if mod(short,2) %short is 1
        xdiff = [xdiff, xlong(closestIndex) - xshort(a)]; % x2-x1
        ydiff = [ydiff, ylong(closestIndex) - yshort(a)]; % x2-x1

    else %short is 2
        xdiff = [xdiff, -xlong(closestIndex) + xshort(a)]; % x1-x2
        ydiff = [ydiff, -ylong(closestIndex) + yshort(a)]; % x1-x2

    end

    if x1cp(a) %channel 1 has stepped, so what was the interhead separation and lifetime
        
    end

    if x2cp(a) %channel 2 has stepped
        
    end

end

figure()

% On-axis Inter-head separation
subplot(1,2,1)
histogram(xdiff)
title('On-axis Separation')

% Off-axis Inter-head separation
subplot(1,2,2)
histogram(ydiff)
title('Off-axis Separation')

% Step size


% Dwell



steps = [];
dwells = [];

fprintf("Now leaving the dichromatic regime \n")




function B = mean_decimate_array(A, decimation_factor)

% Compute the number of complete groups of decimation_factor
numGroups = floor(length(A) / decimation_factor);

% Truncate the array to only include full groups of decimation_factor
truncatedData = A(1:(numGroups * decimation_factor));

% Reshape the truncated array into 5 rows (each column corresponds to a group of decimation_factor points)
reshapedData = reshape(truncatedData, decimation_factor, []);

% Compute the mean of each column
B = mean(reshapedData, 1, 'omitnan');