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
    
for i=1:fnum
    
    if isfile(directory)
        steptrace = load(fname);
    else
        fname = f(i).name;
        steptrace = load(strcat(cd,'\',fname));
    end
    
    trace = steptrace.data;
    if ~isfield(trace,'trace') || ~isfield(trace,'trace_yx')
        fnum = fnum - 1;
    else
    data = trace.trace;
    data_yx = trace.trace_yx;
    
    % Steps
    [on_steps, ~] = add_to_list_6col_steps_v2(data,threshold);
    ONsteps = [ONsteps; on_steps'];
    
    [off_steps, ~] = add_to_list_6col_steps_v2(data_yx,threshold);
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
    
    end %if trace exists
    
end


end

