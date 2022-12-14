% StepandDwellHist.m                    Created by John Canty 03-06-2017

% Edited by Mohamed Elshenawy 07-31-2017
% Then by Joseph Slivka 2022/11/17

% Generates histograms from the steps, dwells, forward and backward dwells from FionaViewer.m

% SubFunctions:
%   add_to_list_6col_steps_v2.m
%   add_to_list_6col_dwells_v2.m
%   add_to_list_6col_dwells_for_back.m;

% Usage:
%   Navigate to the directory containing your .mat files. This function
%   will identify all .mat files in the directory, compile them, then
%   generate the step-size and dwell-time histograms.

function [ONsteps,OFFsteps,dwells,dwells_for,dwells_back] = StepandDwellHist_v2(directory,threshold,framerate)
% Default threshold is 0.

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
    if ~isfield(trace,'trace')
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

% JS Edit 2022/11/11
% Let's pass this into a function. This way we can make it however we want
% and for neighbors also.
% Also option for common titles eventually?

if ~isfile(directory) %only plot for summary
PlotStepStats(fnum, ONsteps, OFFsteps, dwells, dwells_for, dwells_back, fullfile(directory,"AllStats.fig"))

%% If neighbors are a thing we want to examine, run to see neighbor statistics

%Options to put in regions here, maybe if a GUI comes
%xb = [0,50,100,100]; yb = [0,25,25,50];
%xb = [0,100,200]; yb = [0,200,200];
%xb = [0,75,150,225]; yb = [0,200,200,200];
%xb = [0,50,100,200]; yb = [0,200,200,200];
%xa = 0.2*xb; ya = yb;
%Forward/Backward Scheme
xb = [0,50,50,100,100]; yb = [0,25,25,50,50];
xa = [0,0,50,50,100]; ya = yb;

% Generate an automatic foldername that carries Neighbor Info
totarr = [xb,yb,xa,ya];
foldername = '[';
for j = 1:length(totarr)
    foldername = strcat(foldername, num2str(totarr(j)), ',');
    if mod(j,length(xb)) == 0 && j ~= length(totarr)
        foldername = strcat(foldername(1:end-1), '],[');
    end
end
foldername = strcat(foldername(1:end-1), ']');
foldername = fullfile(directory,foldername);
if ~isfolder(foldername)
    mkdir(foldername)
end

fNeighborlyRegions(framerate,directory,xb,yb,xa,ya,0,foldername)

% Finally, compile all the data into some summaries for an excel sheet
StepInfoUnique(framerate,fullfile(directory,'/'),xb,yb,xa,ya,foldername)
% This is a bit recursive, but it never enters this part of the function
% once we set StepInfo in motion

end

