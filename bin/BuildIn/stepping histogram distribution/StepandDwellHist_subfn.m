function [fnum,ONsteps,OFFsteps,dwells,dwells_for,dwells_back] = StepandDwellHist_subfn(directory,threshold,framerate,options)
% Function for StepandDwellHist which goes through all files in a directory
% or just one file if passed and returns the statistics

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
    
    trace = steptrace.data;
    if ~isfield(trace,'trace') %|| ~isfield(trace,'trace_yx')
        fnum = fnum - 1;
    else
    data = trace.trace;
    data_yx = trace.trace; % data_yx = trace.trace_yx;
    %JS Edit 2024/03/07 for loading MINFLUX times rather than framerate
    if isfield(trace,'time')
        framerate = trace.time; %framerate is now actually an array of times
    end

    
    % going to add in option to ignore steps that have too little data
    % points (i.e. <7 to increase accuracy);
        
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
        [forward,backward] = add_to_list_6col_dwells_for_back(data,framerate);
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

