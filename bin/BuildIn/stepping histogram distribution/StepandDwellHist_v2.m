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

function [ONsteps,OFFsteps,dwells,dwells_for,dwells_back] = StepandDwellHist_v2(directory,threshold,framerate,options)
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
<<<<<<< HEAD
    end %if trace exists

    % for pause events
    cnt_pause_events = 47; %make sure to also change in Neighbor Regions
    % Note we are always truncating off the last point in Neighbor Regions
    % since looking for transitions (which require two data points).
    % Therefore we will also truncate exactly one data point per molecule
    % to be consistent
    [pf, Values] = fPauseAnalysis(data(1:end-1,:), [], cnt_pause_events);
    HistVals = [HistVals, Values];
    pause_frequency = [pause_frequency, pf];

=======
    
    end %if trace exists
    
>>>>>>> 2effac3 (GUI integration of Pause and Stepping Analysis with regions)
end

% JS Edit 2022/11/11
% Let's pass this into a function. This way we can make it however we want
% and for neighbors also.
% Also option for common titles eventually?

if ~isfile(directory) %only plot for summary
PlotStepStats(fnum, ONsteps, OFFsteps, dwells, dwells_for, dwells_back, fullfile(directory,"AllStats.fig"))

%% If neighbors are a thing we want to examine, run to see neighbor statistics
if options.UseNeighborRegions
%Options to put in regions here if GUI is not really your thing, but just
%override
%xb = [0,100,200]; yb = [0,200,200];
%xb = [0,75,150,225]; yb = [0,200,200,200];
%xb = [0,50,100,200]; yb = [0,200,200,200];
% xa = xb; ya = yb;
%Forward/Backward Scheme
% xb = 1.5*[0,50,50,100,100]; yb = [0,35,35,70,70];
% xa = 1.5*[0,0,50,50,100]; ya = yb;
<<<<<<< HEAD
=======
% 
% options.XB = xb; options.XA = xa; options.YB = yb; options.YA = ya;
>>>>>>> 2effac3 (GUI integration of Pause and Stepping Analysis with regions)

% Generate an automatic foldername that carries Neighbor Info
totarr = [options.XB,options.YB,options.XA,options.YA];
foldername = '[';
for j = 1:length(totarr)
    foldername = strcat(foldername, num2str(totarr(j)), ',');
    if mod(j,length(options.XB)) == 0 && j ~= length(totarr)
        foldername = strcat(foldername(1:end-1), '],[');
    end
end
foldername = strcat(foldername(1:end-1), ']');
foldername = fullfile(directory,foldername);
if ~isfolder(foldername)
    mkdir(foldername)
end

options.mode = 'steps';
fNeighborlyRegions(options,framerate,directory,0,foldername)


% Finally, compile all the data into some summaries for an excel sheet
StepInfoUnique(options,framerate,fullfile(directory,'/'),foldername)
% This is a bit recursive, but it never enters this part of the function
% once we set StepInfo in motion
<<<<<<< HEAD

% To Set Pause Threshold
figure()
hh = histogram(HistVals,'BinWidth',1);
CSum = zeros(1,hh.NumBins);
CSum(1) = hh.Values(1);
for b = 2:hh.NumBins
    CSum(b) = CSum(b-1) + hh.Values(b);
end
Cnorm = CSum/CSum(end);
% find tau by finding where population is at 1/e
[~,tau] = min(abs(Cnorm - exp(-1)));
% find actual 95% confidence interval by finding where population is at 95%
[~,conf2] = min(abs(Cnorm - 0.99))
hold on
plot((conf2+1)*ones(1,2),[0,max(hh.Values)],'r--');
=======
end

end
>>>>>>> 2effac3 (GUI integration of Pause and Stepping Analysis with regions)

end
