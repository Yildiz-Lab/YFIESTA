function xy_deltatxy_step = two_color_statistics(ch1data,ch2data)
% Take two color data and get statistics:
%   interhead separation and
%   relative rates of stepping
%   step size on interhead separation
%   

decimation_factor = 1;

fprintf("Now entering the dichromatic regime \n")

% Process the data
if isfield(ch1data, 'time')
    time1 = ch1data.time;
    time1 = ch1data.time(~isnan(ch1data.time));
else
    time1 = 1:size(ch1data.xy,1);
end
idx1 = find(~isnan(ch1data.trace(:,1)));

if isfield(ch2data, 'time')
    time2 = ch2data.time;
    time2 = ch2data.time(~isnan(ch2data.time));
else
    time2 = 1:size(ch2data.xy,1);
end
idx2 = find(~isnan(ch2data.trace(:,1)));

trace_start_time = min([time1(1), time2(1)]);

time1 = mean_decimate_array(time1, decimation_factor);
time2 = mean_decimate_array(time2, decimation_factor);

% Only include points in the intersection after we have removed all NaNs
[t_begin, idx_begin] = max([time1(1), time2(1)]); [t_end, idx_end] = min([time1(end), time2(end)]);
if t_begin < 2 
    [~,closestIndex] = min(abs(time2-time1(1)));
    begin_idx_arr = [1, closestIndex];
else
    [~,closestIndex] = min(abs(time1-time2(1)));
    begin_idx_arr = [closestIndex, 1];
end

if t_end < 2 
    [~,closestIndex] = min(abs(time2-time1(end)));
    end_idx_arr = [length(time1), closestIndex];
else
    [~,closestIndex] = min(abs(time1-time2(end)));
    end_idx_arr = [closestIndex, length(time2)];
end

time1 = time1(begin_idx_arr(1):end_idx_arr(1));
time2 = time2(begin_idx_arr(2):end_idx_arr(2));

% Ch1 data %x=on-axis  %y=off-axis
x1 = mean_decimate_array(ch1data.trace(idx1,3), decimation_factor);
y1 = mean_decimate_array(ch1data.trace_yx(idx1,3), decimation_factor);
% And cut off at the appropriate time points
x1 = x1(begin_idx_arr(1):end_idx_arr(1)); y1 = y1(begin_idx_arr(1):end_idx_arr(1));
y1 = -y1; %Need to switch to negative due to wrong cross sign for some reason

% Ch1 data changepoints adapted
x1cp = mean_decimate_array(ch1data.trace(idx1,5), decimation_factor); x1cp = ceil(x1cp);
y1cp = mean_decimate_array(ch1data.trace_yx(idx1,5), decimation_factor); y1cp = ceil(y1cp);
x1cp = x1cp(begin_idx_arr(1):end_idx_arr(1)); y1cp = y1cp(begin_idx_arr(1):end_idx_arr(1));

% Ch2 data  %x=on-axis  %y=off-axis
x2 = mean_decimate_array(ch2data.trace(idx2,3), decimation_factor);
y2 = mean_decimate_array(ch2data.trace_yx(idx2,3), decimation_factor);
x2 = x2(begin_idx_arr(2):end_idx_arr(2)); y2 = y2(begin_idx_arr(2):end_idx_arr(2));

% Ch2 data changepoints adapted
x2cp = mean_decimate_array(ch2data.trace(idx2,5), decimation_factor); x2cp = ceil(x2cp);
y2cp = mean_decimate_array(ch2data.trace_yx(idx2,5), decimation_factor); y2cp = ceil(y2cp);
x2cp = x2cp(begin_idx_arr(2):end_idx_arr(2)); y2cp = y2cp(begin_idx_arr(2):end_idx_arr(2));

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
end

% Now to do time, we should make a time organized list with all of the
% information in rows. Then, we just walk through.

% The matrix (T_org) should be the following:
% Col 1: Time (tim) % organized in ascending order
% Col 2: x ch1 (on-axis)
% Col 3: y ch1 (off-axis)
% Col 4: x ch2
% Col 5: y ch2
% Col 6: xcp ch1 (on-axis changepoint)
% Col 7: ycp ch1 (off-axis changepoint)
% Col 8: xcp ch2 (on-axis changepoint)
% Col 9: ycp ch2 (on-axis changepoint)

T_org = [time1', x1', y1', nan(length(time1),1), nan(length(time1),1), x1cp', y1cp', zeros(length(time1),1), zeros(length(time1),1)];
T_org = [T_org; time2', nan(length(time2),1), nan(length(time2),1), x2', y2', zeros(length(time2),1), zeros(length(time2),1), x2cp', y2cp'];
T_org = sortrows(T_org, 1);

% [(1:200)', T_org(1:200,1:5)]

ch1mask = find(~isnan(T_org(:,2)));
ch2mask = find(~isnan(T_org(:,4)));

ch1xcpmask = find(T_org(:,6));
ch1ycpmask = find(T_org(:,7));
ch2xcpmask = find(T_org(:,8));
ch2ycpmask = find(T_org(:,9));

xy_deltatxy_step = zeros(0,7); %store: channel, time, sep x, sep y, step delta t, step delta x, step delta y
% debugging_matrix = zeros(0,7);
last_ch1_idx = 1;
last_ch2_idx = 1;

for a = 1:size(T_org,1) % could actually probably just do intersect changepoints for this one if we wanted to...
    if T_org(a,6) %ch1 x-changepoint, ch1 x-step
        
        am1 = max(intersect(ch1mask,1:a-1)); % to get nearest point one behind in the same channel
        
        [~,closestIndex] = min(abs(T_org(ch2mask,1) - T_org(a,1))); % find nearest in ch2 (honestly, we could just search the neighborhood. Later)
        closestIndex = ch2mask(closestIndex); %convert to the actual index rather than mask index

        [~,closestIndex_m1] = min(abs(T_org(ch2mask,1) - T_org(am1,1))); % find nearest in ch2 one point behind
        closestIndex_m1 = ch2mask(closestIndex_m1); %convert to the actual index rather than mask index

        % define as conditioned on the state of the other head (i.e. final
        % - initial)
        stepxdiff = mean(T_org(last_ch2_idx:closestIndex_m1,4),'omitnan') - mean(T_org(last_ch1_idx:am1,2), 'omitnan');
        stepydiff = mean(T_org(last_ch2_idx:closestIndex_m1,5),'omitnan') - mean(T_org(last_ch1_idx:am1,3), 'omitnan');
        
        ap1 = min(intersect(ch1xcpmask,a+1:size(T_org,1))); % to get next nearest changepoint ahead in the same channel
        ap1 = max(intersect(ch1mask,a:ap1-1));

        deltat = T_org(a,1) - T_org(last_ch1_idx,1);
        deltax = mean(T_org(a:ap1,2), 'omitnan') - mean(T_org(last_ch1_idx:am1,2), 'omitnan');
        deltay = NaN; %mean(T_org(a:ap1,3), 'omitnan') - mean(T_org(last_ch1_idx:am1,3), 'omitnan');
        
        % debugging_matrix = [debugging_matrix; 1, last_ch1_idx, am1, a, last_ch2_idx, closestIndex_m1, closestIndex];

        % Update storage
        new_row = [1, T_org(a,1)-trace_start_time, stepxdiff, stepydiff, deltat, deltax, deltay];
        xy_deltatxy_step = [xy_deltatxy_step; new_row];
        % and reset the clock
        last_ch1_idx = a;
        last_ch2_idx = closestIndex;
    end

    if T_org(a,7) %ch1 y-changepoint, ch1 y-step

        am1 = max(intersect(ch1mask,1:a-1)); % to get nearest point one behind in the same channel

        [~,closestIndex] = min(abs(T_org(ch2mask,1) - T_org(a,1))); % find nearest in ch2 (honestly, we could just search the neighborhood. Later)
        closestIndex = ch2mask(closestIndex); %convert to the actual index rather than mask index

        [~,closestIndex_m1] = min(abs(T_org(ch2mask,1) - T_org(am1,1))); % find nearest in ch2 one point behind
        closestIndex_m1 = ch2mask(closestIndex_m1); %convert to the actual index rather than mask index

        % define as conditioned on the state of the other head (i.e. final
        % - initial)
        stepxdiff = mean(T_org(last_ch2_idx:closestIndex_m1,4),'omitnan') - mean(T_org(last_ch1_idx:am1,2), 'omitnan');
        stepydiff = mean(T_org(last_ch2_idx:closestIndex_m1,5),'omitnan') - mean(T_org(last_ch1_idx:am1,3), 'omitnan');

        ap1 = min(intersect(ch1ycpmask,a+1:size(T_org,1))); % to get next nearest changepoint ahead in the same channel
        ap1 = max(intersect(ch1mask,a:ap1-1));

        deltat = T_org(a,1) - T_org(last_ch1_idx,1);
        deltax = NaN; %mean(T_org(a:ap1,2), 'omitnan') - mean(T_org(last_ch1_idx:am1,2), 'omitnan');
        deltay = mean(T_org(a:ap1,3), 'omitnan') - mean(T_org(last_ch1_idx:am1,3), 'omitnan');

        % debugging_matrix = [debugging_matrix; 1, last_ch1_idx, am1, a, last_ch2_idx, closestIndex_m1, closestIndex];

        % Update storage
        new_row = [1, T_org(a,1)-trace_start_time, stepxdiff, stepydiff, deltat, deltax, deltay];
        xy_deltatxy_step = [xy_deltatxy_step; new_row];
        % and reset the clock
        last_ch1_idx = a;
        last_ch2_idx = closestIndex;
    end

    if T_org(a,8) %ch2 x-changepoint, ch2 x-step
        
        am1 = max(intersect(ch2mask,1:a-1)); % to get nearest point one behind in the same channel
        
        [~,closestIndex] = min(abs(T_org(ch1mask,1) - T_org(a,1))); % find nearest in ch2 (honestly, we could just search the neighborhood. Later)
        closestIndex = ch1mask(closestIndex); %convert to the actual index rather than mask index
        
        [~,closestIndex_m1] = min(abs(T_org(ch1mask,1) - T_org(am1,1))); % find nearest in ch2 one point behind
        closestIndex_m1 = ch1mask(closestIndex_m1); %convert to the actual index rather than mask index

        % define as conditioned on the state of the other head (i.e. final
        % - initial)
        stepxdiff = mean(T_org(last_ch1_idx:closestIndex_m1,2),'omitnan') - mean(T_org(last_ch2_idx:am1,4), 'omitnan');
        stepydiff = mean(T_org(last_ch1_idx:closestIndex_m1,3),'omitnan') - mean(T_org(last_ch2_idx:am1,5), 'omitnan');
        
        ap1 = min(intersect(ch2xcpmask,a+1:size(T_org,1))); % to get next nearest changepoint ahead in the same channel
        ap1 = max(intersect(ch2mask,a:ap1-1));

        deltat = T_org(a,1) - T_org(last_ch2_idx,1);
        deltax = mean(T_org(a:ap1,4), 'omitnan') - mean(T_org(last_ch2_idx:am1,4), 'omitnan');
        deltay = NaN; %mean(T_org(a:ap1,5), 'omitnan') - mean(T_org(last_ch2_idx:am1,5), 'omitnan');
        
        % debugging_matrix = [debugging_matrix; 2, last_ch2_idx, am1, a, last_ch1_idx, closestIndex_m1, closestIndex];

        % Update storage
        new_row = [2, T_org(a,1)-trace_start_time, stepxdiff, stepydiff, deltat, deltax, deltay];
        xy_deltatxy_step = [xy_deltatxy_step; new_row];
        % and reset the clock
        last_ch1_idx = closestIndex;
        last_ch2_idx = a;
    end

    if T_org(a,9) %ch2 y-changepoint, ch2 y-step
        
        am1 = max(intersect(ch2mask,1:a-1)); % to get nearest point one behind in the same channel
        
        [~,closestIndex] = min(abs(T_org(ch1mask,1) - T_org(a,1))); % find nearest in ch2 (honestly, we could just search the neighborhood. Later)
        closestIndex = ch1mask(closestIndex); %convert to the actual index rather than mask index
        
        [~,closestIndex_m1] = min(abs(T_org(ch1mask,1) - T_org(am1,1))); % find nearest in ch2 one point behind
        closestIndex_m1 = ch1mask(closestIndex_m1); %convert to the actual index rather than mask index

        % define as conditioned on the state of the other head (i.e. final
        % - initial)
        stepxdiff = mean(T_org(last_ch1_idx:closestIndex_m1,2),'omitnan') - mean(T_org(last_ch2_idx:am1,4), 'omitnan');
        stepydiff = mean(T_org(last_ch1_idx:closestIndex_m1,3),'omitnan') - mean(T_org(last_ch2_idx:am1,5), 'omitnan');
        
        ap1 = min(intersect(ch2ycpmask,a+1:size(T_org,1))); % to get next nearest changepoint ahead in the same channel
        ap1 = max(intersect(ch2mask,a:ap1-1));

        deltat = T_org(a,1) - T_org(last_ch2_idx,1);
        deltax = NaN; %mean(T_org(a:ap1,4), 'omitnan') - mean(T_org(last_ch2_idx:am1,4), 'omitnan');
        deltay = mean(T_org(a:ap1,5), 'omitnan') - mean(T_org(last_ch2_idx:am1,5), 'omitnan');
        
        % debugging_matrix = [debugging_matrix; 2, last_ch2_idx, am1, a, last_ch1_idx, closestIndex_m1, closestIndex];

        % Update storage
        new_row = [2, T_org(a,1)-trace_start_time, stepxdiff, stepydiff, deltat, deltax, deltay];
        xy_deltatxy_step = [xy_deltatxy_step; new_row];
        % and reset the clock
        last_ch1_idx = closestIndex;
        last_ch2_idx = a;
    end

end

% xdiff
% ydiff
% xy_deltatxy_step
% debugging_matrix

figure()
% On-axis Inter-head separation
subplot(1,2,1)
histogram(xdiff)
title('On-axis Separation (nm)')

% Off-axis Inter-head separation
subplot(1,2,2)
histogram(ydiff)
title('Off-axis Separation (nm)')

figure()
% Step size
subplot(3,1,1)
scatter(xy_deltatxy_step(:,3), xy_deltatxy_step(:,6), 50, 'k', 'filled')
ylabel('\Delta x')
title('Dependence on inter-head separation')

subplot(3,1,2)
scatter(xy_deltatxy_step(:,4), xy_deltatxy_step(:,7), 50, 'k', 'filled')
ylabel('\Delta x')
title('Dependence on inter-head separation')

% Dwell
subplot(3,1,3)
scatter(xy_deltatxy_step(:,3), xy_deltatxy_step(:,5), 50, 'k', 'filled')
ylabel('Dwell time')
xlabel('Head Separation (nm)')

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