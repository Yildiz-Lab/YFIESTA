% Plot all stats from a selection of MATLAB traces as figures.
% Options add the ability to make a png or pdf plot

% BEWARE CHANNELS GET EASILY MIXED UP!! MAKE SURE TO LABEL CORRECTLY CH1 OR
% CH2 WITH YOUR MEASUREMENT BEFORE SAVING!!!

% For purposes of plotting, Ch1 is the first channel and Ch2 is the
% following neighbor
% For purposes of coding, all subscript 1's go with Ch2 and all subscript
% 2's go with Ch1. I mixed them up at some point.

offset = -15;
dx = -18; %dx = 75; %assume on-axis, horizontal MT
dy = -27; %dy = 35; %assume off-axis, vertical MT
% Note, can instead use function offset_by_law_of_means(dir, fnames) below

% Options
grid_spacing = 16;
%number of points external to mean/smooth with on both ends (so actually doubled).
% Thus, a sraw = 2 implies you are using 5 points, the original point and two on either side.
% For no running average change set to 0.
smooth_running_average_wdw = 0;


[fnames, mydirectory] = uigetfile({'*.mat'},'MultiSelect','on');
%Multiselect allows user to grab more than one file and parse accordingly

% if only one file is selected, need to package into a cell array so that
% the for loop will work
if class(fnames) == 'char'
    fnames = {fnames};
end

% option to set offsets for plots!
% [dx, dy] = offset_by_law_of_means(dir, fnames);
dx = 0;
dy = 0;
offset = -dx+40; %arbitrary value so that on and off axes don't overlap

spatiotemporal_info = zeros(2,2,1);

% loop through each selected file, assuming they are all in the same
% directory
    
for i = 1:length(fnames)
    try
        fname = fnames{i};
    catch
        fprintf("Finished combining xy and yx traces. Aborting... \n")
        return
    end

    % MOLECULE
    % find the molecule in the data Molecule if loaded
    jj = strfind(fname,'_fiona')-1; %nbh or initial
    j = strfind(fname,'_nbh')-1; %if there is a neighbor
    ii = 1; gg = 1;
    if ~isempty(j)
        ii = j+6;
    end

    if contains(fname(1:j),'_') % actually a more complicated name
        gg = strfind(fname(1:j),'_')+1;
        if length(gg) > 1
        gg = gg(end);
        end
    end
    
    % fname(ii:jj)
    % fname(gg:j)

    if ~isempty(j) %switch it

        fname2 = fname;

        strstart = strfind(fname,'nbh');
        strend = strfind(fname,'fiona');
        fname1 = strcat(fname(1:strstart-1),fname(strend:end));
    else
        fname1 = fname;
    end
    
    fname1
    fname2

    % Now dissect filename connections
    totdata = load(fullfile(mydirectory,fname1));
    ch1data = totdata.data;

    if isfield(ch1data, 'time')
        time1 = ch1data.time;
        time1 = ch1data.time(~isnan(ch1data.time));
    else
        time1 = 1:size(ch1data.xy,1);
    end

    idx1 = find(~isnan(ch1data.trace(:,1)));

    % Smooth by mean window
    ch1trace = smooth_running_average(ch1data.trace(idx1,1), smooth_running_average_wdw);
    ch1trace_yx = smooth_running_average(ch1data.trace_yx(idx1,1), smooth_running_average_wdw);


    if ~isempty(j) % now do the neighbor

        % try

            totdata = load(fullfile(mydirectory,fname2));
            ch2data = totdata.data;
    
            if isfield(ch2data, 'time')
                time2 = ch2data.time;
                time2 = ch2data.time(~isnan(ch2data.time));
            else
                time2 = 1:size(ch2data.xy,1);
            end
    
            idx2 = find(~isnan(ch2data.trace(:,1)));
            % time1 is not a typo, want to make sure you are normalizing to the
            % same time

            % Smooth by mean window
            ch2trace = smooth_running_average(ch2data.trace(idx2,1), smooth_running_average_wdw);
            ch2trace_yx = smooth_running_average(ch2data.trace_yx(idx2,1), smooth_running_average_wdw);

    else
        fprintf("no associated neighbor file ... continuing with just a single channel \n")

        idx2 = [NaN];

    end
    
    st = zeros(2,2);

    if sum(~isnan(idx2)) > 0

        t1 = max([time1(1),time2(1)]);
        t2 = min([time1(end),time2(end)]);

        midx1a = find(time1 < t2); midx1b = find(time1 > t1);
        midx2a = find(time2 < t2); midx2b = find(time2 > t1);
        idx1 = idx1(intersect(midx1a,midx1b));
        idx2 = idx2(intersect(midx2a,midx2b));

        y1 = (ch1data.trace_yx(idx1(2:end),1) - ch1data.trace_yx(idx1(1:end-1),1)).^2;
        x1 = (ch1data.trace(idx1(2:end),1) - ch1data.trace(idx1(1:end-1),1) ).^2;

        st(1,1) = mean( sqrt(x1 + y1)) / sqrt(2);
        st(1,2) = mean(time1(2:end)-time1(1:end-1))*1000;

        fprintf(strcat("ch1 dx: ", num2str( round(st(1,1),2) ), " nm \n"))
        fprintf(strcat("ch1 dt: ", num2str( round(st(1,2),3) ), " ms \n"))

        y2 = (ch2data.trace_yx(idx2(2:end),1) - ch2data.trace_yx(idx2(1:end-1),1)).^2;
        x2 = (ch2data.trace(idx2(2:end),1) - ch2data.trace(idx2(1:end-1),1) ).^2;
        
        st(2,1) = mean( sqrt(x2 + y2)) / sqrt(2);
        st(2,2) = mean(time2(2:end)-time2(1:end-1))*1000;
    
        fprintf(strcat("ch2 dx: ", num2str( round(st(2,1),2)), " nm \n"))
        fprintf(strcat("ch2 deltat: ", num2str( round(st(2,2),3) ), " ms \n"))
    
    else

        y1 = (ch1data.trace_yx(idx1(2:end),1) - ch1data.trace_yx(idx1(1:end-1),1)).^2;
        x1 = (ch1data.trace(idx1(2:end),1) - ch1data.trace(idx1(1:end-1),1) ).^2;
        
        st(1,1) = mean( sqrt(x1 + y1)) / sqrt(2);
        st(1,2) = mean(time1(2:end)-time1(1:end-1))*1000;
    
        fprintf(strcat("ch1 dx: ", num2str( round(st(1,1),2)), " nm \n"))
        fprintf(strcat("ch1 deltat: ", num2str( round(st(1,2),3) ), " ms \n"))


    end

    spatiotemporal_info(:,:,i) = st;

end

spatiotemporal_info
average_spatiotemporal_info = mean(spatiotemporal_info,3)

figure()
subplot(2,2,1)
histogram(spatiotemporal_info(1,1,:))
title('Channel 1 \Delta x')
legend(strcat('\mu = ', num2str(round(average_spatiotemporal_info(1,1),2))))
subplot(2,2,2)
histogram(spatiotemporal_info(1,2,:))
title('Channel 1 \Delta t')
legend(strcat('\mu = ', num2str(round(average_spatiotemporal_info(1,2),2))))
subplot(2,2,3)
histogram(spatiotemporal_info(2,1,:))
title('Channel 2 \Delta x')
legend(strcat('\mu = ', num2str(round(average_spatiotemporal_info(2,1),2))))
subplot(2,2,4)
histogram(spatiotemporal_info(2,2,:))
title('Channel 2 \Delta t')
legend(strcat('\mu = ', num2str(round(average_spatiotemporal_info(2,2),2))))


function A = smooth_running_average(A, smooth_running_average_wdw)

if smooth_running_average_wdw < 1
    return
else

ll = length(A);
A = [nan(smooth_running_average_wdw,1); A; nan(smooth_running_average_wdw,1)];
for c = 1:ll
    A(c) = mean(A(c:c+2*smooth_running_average_wdw),'omitnan');
end
% now delete the excess
A(end+1-2*smooth_running_average_wdw:end) = [];

end

end