function [fnum,ONsteps,OFFsteps,dwells,dwells_for,dwells_back] = StepandDwellHist_subfn(directory,threshold,framerate,options)
% Function for StepandDwellHist which goes through all files in a directory
% or just one file if passed and returns the statistics

framerate = 0.1;
% check if directory is a file or a folder
if isfile(directory)
    fname = directory;
    fnum = 1;
else
    % Gather steps and dwells all in one folder
    cd = directory; %JS Edit 220207
    f = dir(fullfile(cd,'*.mat')); %JS Edit 220207
    fnum = length(f);
end

% Make blank arrays for storage of 1C ON, OFFsteps and dwells forward and
% backward
ONsteps = [];
OFFsteps = [];
dwells = [];
dwells_for = [];
dwells_back = [];

% Also blank arrays for polar plots
r = [];
theta = [];

% Blank multipurpose array for two_color_statistics
two_color_xy_deltatxy_step = zeros(0,7);
two_color_xydiff = zeros(0,2);

if nargin < 2
    disp("Forgotten threshold or framerate")
    return
end

% Don't know if this is still necessary for MAC but doesn't hurt
if fnum > 1 
    for i=fnum:-1:1
        if contains(f(i).name,'._') %JS Edit to ignore extra '._' that randomly show up sometimes MAC
        f(i) = [];
        end
    end
    fnum = length(f);
else
    if contains(fname,'._') %JS Edit to ignore extra '._' that randomly show up sometimes MAC
        fnum = 0;
        return
    end
end

% Loop through each file in the directory
for i=1:fnum
    
    if isfile(directory) %if the directory is actually a file (i.e. only one)
        [d,f,e] = fileparts(directory);
        % if contains(f,'._') %JS Edit to delete extra ._ from an error
        % f = f(3:end);
        % end
        fname = fullfile(d,strcat(f,e));
        steptrace = load(fullfile(d,strcat(f,e)));
    else
        fname = f(i).name;
        % if contains(fname,'._') %JS Edit to delete extra ._ from an error
        % fname = fname(3:end);
        % end
        steptrace = load(fullfile(cd,'/',fname));
    end
    
    %% Processing Trace into steps
    % Extract trace from data structure
    if isfield(steptrace,'data')
        trace = steptrace.data;
    elseif isfield(steptrace,'trace')
        trace = steptrace;
    end
    
    % Skip if trace does not exist
    if ~isfield(trace,'trace') || ~isfield(trace,'trace_yx')
        fnum = fnum - 1;
    else

    %JS Edit 2024/03/07 for loading MINFLUX times rather than framerate
    if isfield(trace,'time')
        framerate = trace.time; %framerate is now actually an array of times
        framerate = framerate(~isnan(framerate)); % also ignoring Nans
    end
    
    % % OLD FILTERING (NOW NICELY PUT INTO FUNCTIONS)
    % data = trace.trace;
    % data_yx = trace.trace_yx;
    % 
    % % going to add in option to ignore steps that have too little data
    % % points (i.e. <7 to increase accuracy);
    % 
    % % actually adjust to have it not just remove the steps but add it to
    % % the step that leads to the smallest change in position, maybe because
    % % they are just short backwards dips
    % remove_dips_time = str2double(options.OmitBlips)/1000; % convert to seconds
    % data = modify_filter_trace(data,remove_dips_time,[],framerate);
    % data_yx = modify_filter_trace(data_yx,remove_dips_time,[],framerate);
    % trace.trace = data; trace.trace_yx = data_yx;
    % 
    % % JS Edit 2025/01/08
    % % Get changepoints back in for 2C data with many Nans.
    % % Should have been careful changing the changepoint code.
    % nnidx = find(~isnan(data(:,1)));
    % 
    % % on-axis
    % data = data(nnidx,:);
    % chp = find( abs(data(2:end,3) - data(1:end-1,3)) > 0);
    % data(chp+1,5) = 1;
    % 
    % % off-axis
    % data_yx = data_yx(nnidx,:);
    % chp = find( abs(data_yx(2:end,3) - data_yx(1:end-1,3)) > 0);
    % data_yx(chp+1,5) = 1;
    % 
    % 
    % % if want to merge_step_components in post
    % if options.Merge
    %     % automatically show polar plots
    %     [rprime, thetaprime] = polar_conversion(trace,0);
    %     r = [r, rprime]; theta = [theta, thetaprime];
    % 
    %     trace = merge_step_components(trace);
    %     data = trace.trace_2d;
    %     data_yx = trace.trace_2d;
    % end
    
    % JS Edit 2025/02/27 just making it easier to filter in the future.
    % More compartmentalized

    % Convert trace
    remove_dips_time = str2double(options.OmitBlips)/1000; % convert to seconds
    trace = filterSteps(trace, framerate, remove_dips_time);
    
    % Merge if desired, also return r, thetha
    if options.Merge
        trace = mergeSteps(trace);
        r = [r, trace.r];
        theta = [theta, trace.theta];
        data = trace.trace_2d;
        data_yx = trace.trace_2d;
    else
        data = trace.trace;
        data_yx = trace.trace_yx;
    end
    
    %% Extract data from prepared steps
    % Steps
    [on_steps, ~] = add_to_list_6col_steps_v2(data,threshold);
    ONsteps = [ONsteps; on_steps'];
    
    %[off_steps, ~] = add_to_list_6col_steps_v2(data_yx,threshold); %using complete mean of step
    [off_steps, ~] = add_to_list_6col_steps_v2_instant(data_yx,threshold); %using a mean of a smaller window (a derivative). Change window size in function.
    OFFsteps = [OFFsteps; off_steps'];
    
    % Dwells
    mat = add_to_list_6col_dwells_v2(data,threshold,[],framerate);
    if ~isempty(mat) %check that dwells were found (JS Edit 220310)
        dwell = mat(:,3);
        dwells = [dwells; dwell];

        % Forward and Backward dwells
        if options.TieDwells % tie to the previous step (so the dwell is after the step)
        [forward,backward] = add_to_list_6col_dwells_after_for_back(data,framerate);
        else % otherwise normally where the ending is the rate, this is the common and actual definition way to do it
            [forward,backward] = add_to_list_6col_dwells_for_back(data,framerate);
        end
        dwells_for = [dwells_for; forward];
        dwells_back = [dwells_back; backward];
    end

    
    % JS Edit 2025/02/27  2-Color Statistics
    % if the file is a neighbor file, then we will look in the same
    % directory for its associated data file. Then we will send them both
    % to extract the 2C stepping statistics

    % MOLECULE
    % find the molecule in the data Molecule if loaded
    jj = strfind(fname,'_fiona')-1; %nbh or initial
    j = strfind(fname,'_nbh')-1; %if there is a neighbor
    ii = 1;
    if ~isempty(j)
        ii = j+6;
    end
    
    % This is a neighbor, part of a pair, called a neighbor.
    % We should now look for our other data parent and load both datasets.
    if ~isempty(j) 
        
        % We already did all the filtering for the first stage
        fname2 = fname;
        ch2trace = trace;

        strstart = strfind(fname,'nbh');
        strend = strfind(fname,'fiona');
        % JS Edit works for filename types of
        % (mol#)_nbh_(mol#)_(whatever_whatever)_fiona.mat
        aa = strfind(fname,'_');
        aa = aa(aa > strstart);
        strend = aa(2)+1; %Pull second since first is to prelude the nbh number, then following the neighbor number
        fname1 = strcat(fname(1:strstart-1),fname(strend:end));
    
        % Now dissect filename connections
        totdata = load(fullfile(directory,fname1));
        ch1trace = totdata.data;
        
        ch1trace = filterSteps(ch1trace, framerate, remove_dips_time);
        
        % Now that we have both, let's send the trace data into
        % two_color_statistics
        
        % Also do we want to use merged steps? If so, I need to find out
        % how to pass it
        % Maybe I just change trace before I pass it

        if options.Merge
            ch1trace = mergeSteps(ch1trace);
            ch1trace.trace = ch1trace.trace_2d;
            ch1trace.trace_yx = ch1trace.trace_2d;
            ch2trace.trace = ch2trace.trace_2d;
            ch2trace.trace_yx = ch2trace.trace_2d;
        end
        
        % xy_deltatxy_step (array) - each containing information about a changepoint step with the following format:
        % [channel, normalized trace time (s), x_head_separation (nm), y_head_separation (nm), deltat (s), deltax (nm), deltay (nm)];
        % ch1trace.trace(1:1000,:)
        % ch2trace.trace(1:1000,:)
        two_color_data_struct = two_color_statistics(ch1trace, ch2trace);
        % two_color_data_struct.xydiff(1:1000)
        % two_color_data_struct.xy_deltatxy_step
        two_color_xydiff = [two_color_xydiff; two_color_data_struct.xydiff];
        two_color_xy_deltatxy_step = [two_color_xy_deltatxy_step; two_color_data_struct.xy_deltatxy_step];

    end


%     % JS Edit 2024/03/07
%     %post processing too small dwell steps (incrase accuracy)
%     % find steps below the time threshold. For MINFLUX resolving
%     % ~2ms/point, I will say 0.015s.
%     time_threshold = 0.015;
%     short_steps = find(dwells < time_threshold);
%     
%     % Now that this lot (dwells) should not be used for steps (we can still use it for dwells), we have to remove the fences on either
%     % side which is the i'th and i-1'th steps. Further we should set them
%     % to Nan, and remove after we've gone through
%     
%     for k=length(short_steps):-1:1
%         ONsteps(short_steps(k)) = NaN;
%         if short_steps(k)-1 > 0
%         ONsteps(short_steps(k)-1) = NaN;
%         end
%     end
%     % NaN's will be removed in the final figure
%     % End of JS Edit 2024/03/07
%     
%     end %if trace exists
    
    end


end

if options.Merge && ~isfile(directory)
    plot_polar_conversion(r, theta, 24)
end

% JS Edit 2025/03/05 add in two-color plotting option
if ~isempty(two_color_xydiff)
    Plot2CStepStats(two_color_xydiff, two_color_xy_deltatxy_step)
end


function filtered_trace = filterSteps(trace, framerate, remove_dips_time)
% trace: loaded struct with important fields such as time, xy, yx, trace,
% trace_yx, and neighbors. Will add trace_2d after processing with this
% script
% framerate: a pass in if time is not used. In reality, for changing steps,
% it is not really used
% remove_dips_time: remove small steps below a user defined threshold

    data = trace.trace;
    data_yx = trace.trace_yx;

    % going to add in option to ignore steps that have too little data
    % points (i.e. <7 to increase accuracy);

    % actually adjust to have it not just remove the steps but add it to
    % the step that leads to the smallest change in position, maybe because
    % they are just short backwards dips
    data = modify_filter_trace(data,remove_dips_time,[],framerate);
    data_yx = modify_filter_trace(data_yx,remove_dips_time,[],framerate);

    % JS Edit 2025/01/08
    % Get changepoints back in for 2C data with many Nans.
    % Should have been careful changing the changepoint code.
    nnidx = find(~isnan(data(:,1)));
    
    min_delta_step = 0;
    % on-axis
    data = data(nnidx,:);
    chp = find( abs(data(2:end,3) - data(1:end-1,3)) > min_delta_step);
    data(nnidx,5) = 0; data(chp+1,5) = 1;
    
    % off-axis
    data_yx = data_yx(nnidx,:);
    chp = find( abs(data_yx(2:end,3) - data_yx(1:end-1,3)) > min_delta_step);
    data_yx(nnidx,5) = 0; data_yx(chp+1,5) = 1;

    filtered_trace = trace;
    filtered_trace.trace = data; filtered_trace.trace_yx = data_yx;
    


function merged_trace = mergeSteps(trace)
    % merge_step_components in post
    % automatically show polar plots
    [rprime, thetaprime] = polar_conversion(trace,0); %second is opt to plot invdividual
    % r = [r, rprime]; theta = [theta, thetaprime];

    merged_trace = trace;
    merged_trace = merge_step_components(trace);
    merged_trace.r = rprime;
    merged_trace.theta = thetaprime;