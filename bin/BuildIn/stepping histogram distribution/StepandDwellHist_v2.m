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

options
[fnum,ONsteps,OFFsteps,dwells,dwells_for,dwells_back] = StepandDwellHist_subfn(directory,threshold,framerate,options);

% JS Edit 2022/11/11
% Let's pass this into a function. This way we can make it however we want
% and for neighbors also.
% Also option for common titles eventually?

if ~isfile(directory) %only plot for summary
PlotStepStats(fnum, ONsteps, OFFsteps, dwells, dwells_for, dwells_back, options, fullfile(directory,"AllStats.fig"))
StepInfoUnique(options,framerate,fullfile(directory,'/'))
end

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
% 
% options.XB = xb; options.XA = xa; options.YB = yb; options.YA = ya;

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
end

end

