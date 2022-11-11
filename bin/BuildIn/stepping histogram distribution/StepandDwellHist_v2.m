% StepandDwellHist.m                    Created by John Canty 03-06-2017

% Edited by Mohamed Elshenawy 07-31-2017

% Generates histograms from the steps, dwells, forward and backward dwells from FionaViewer.m

% SubFunctions:
%   add_to_list_6col_steps_v2.m
%   add_to_list_6col_dwells_v2.m
%   add_to_list_6col_dwells_for_back.m;

% Usage:
%   Navigate to the directory containing your .mat files. This function
%   will identify all .mat files in the directory, compile them, then
%   generate the step-size and dwell-time histograms.

function [steps,dwells,dwells_for,dwells_back] = StepandDwellHist_v2(directory,threshold,framerate)
% Default threshold is 0.

% Gather steps and dwells all in one folder
cd = directory; %JS Edit 220207
f = dir(fullfile(cd,'*.mat')); %JS Edit 220207
fnum = length(f);

steps = [];
dwells = [];
dwells_for = [];
dwells_back = [];

if nargin < 2
    disp("Forgotten threshold or framerate")
    return
end
    
for i=1:fnum
    fname = f(i).name;
    steptrace = load(strcat(cd,'\',fname));
    trace = steptrace.data;
    data = trace.trace;
    
    % Steps
    [on_steps, off_steps] = add_to_list_6col_steps_v2(data,threshold);
    steps = [steps;[on_steps', off_steps']];

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
    
end

% JS Edit 2022/11/11
% Let's pass this into a function. This way we can make it however we want
% and for neighbors also.
% Also option for common titles eventually?

PlotStepStats(steps(:,1), steps(:,2), dwells, dwells_for, fnum)



%% If neighbors are a thing we want to examine, run to see neighbor statistics

fNeighborlyRegions(framerate,directory) %Options to put in regions here, maybe if a GUI comes

