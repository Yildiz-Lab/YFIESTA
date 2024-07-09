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

ONsteps = [];
OFFsteps = [];
dwells = [];
dwells_for = [];
dwells_back = [];

if nargin < 2
    disp("Forgotten threshold or framerate")
    return
end

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

r = [];
theta = [];

for i=1:fnum
    
    if isfile(directory)
        [d,f,e] = fileparts(directory);
        % if contains(f,'._') %JS Edit to delete extra ._ from an error
        % f = f(3:end);
        % end
        steptrace = load(fullfile(d,strcat(f,e)));
    else
        fname = f(i).name;
        % if contains(fname,'._') %JS Edit to delete extra ._ from an error
        % fname = fname(3:end);
        % end
        steptrace = load(fullfile(cd,'/',fname));
    end
    
    if isfield(steptrace,'data')
    trace = steptrace.data;
    elseif isfield(steptrace,'trace')
    trace = steptrace;
    end

    if ~isfield(trace,'trace') || ~isfield(trace,'trace_yx')
        fnum = fnum - 1;
    else
    data = trace.trace;
    data_yx = trace.trace_yx;
    
    %JS Edit 2024/03/07 for loading MINFLUX times rather than framerate
    if isfield(trace,'time')
        framerate = trace.time; %framerate is now actually an array of times
    end

    % going to add in option to ignore steps that have too little data
    % points (i.e. <7 to increase accuracy);
    % actually adjust to have it not just remove the steps but add it to
    % the step that leads to the smallest change in position, maybe because
    % they are just short backwards dips
    remove_dips_time = str2double(options.OmitBlips)/1000; % convert to seconds
    data = modify_filter_trace(data,remove_dips_time,[],framerate);
    data_yx = modify_filter_trace(data_yx,remove_dips_time,[],framerate);
    trace.trace = data; trace.trace_yx = data_yx;
    
    % if want to merge_step_components in post
    if options.Merge
        % automatically show polar plots
        [rprime, thetaprime] = polar_conversion(trace,0);
        r = [r, rprime]; theta = [theta, thetaprime];

        trace = merge_step_components(trace);
        data = trace.trace_2d;
        data_yx = trace.trace_2d;
    end

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