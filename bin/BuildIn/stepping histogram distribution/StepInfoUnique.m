function StepInfoUnique(options, framerate, rootfolder, neighborsavefolder)

% JS 2022/11/27
% Compile Individual Collective Data for Traces to help summarize
% Information. Highly versatile and can change to whims of user by just
% changing values in the table.

% All the options to override if desired
if nargin < 3
    rootfolder = uigetdir();
    rootfolder = fullfile(rootfolder,'/');
end

if nargin < 4
    neighborsavefolder = [];
end

f = dir(fullfile(rootfolder,'*.mat')); %JS Edit 220207
fnum = length(f);

% StatsArr
% NumOn, Mean, Std, NumBack, BackMean, BackStd, NumOff, AbsMean, AbsStd
StatsArr = NaN(fnum, 10);
fnames = cell(fnum,1);

for i = 1:fnum
    fnames{i} = f(i).name;
    fullname = strcat(rootfolder,f(i).name);
    StatsArr = FillStatsArray(StatsArr, fullname, i, framerate, options); 
end

StatsTable = ConvertToTable(StatsArr, fnames);
writetable(StatsTable, strcat(rootfolder,'AllStepInfo.xlsx'));

%% NOW TO DO NEIGHBORS IF EXIST
if options.UseNeighborRegions
options.mode = 'steps';
xa = options.XA; xb = options.XB; ya = options.YA; yb = options.YB;

fNeighborlyRegions(options, framerate, rootfolder, 1)
% The Neighbor Data

NStatsArr = NaN(fnum, 10, length(xb));
for i = 1:fnum
    fnames{i} = f(i).name;
    fullname = strcat(rootfolder,f(i).name);
    NStatsArr = FillNeighborStatsArray(NStatsArr, fullname, i, framerate);
end

skipped = 0;
for k = 1:length(xb)
    if any(NStatsArr(:,1,k) > 0)
        if k == length(xb) && skipped == length(xb)-1
            return
        else
            StatsTable = ConvertToTable(NStatsArr(:,:,k), fnames);
            writetable(StatsTable,strcat(rootfolder,'AllNeighborStepInfo.xlsx'),'Sheet',strcat('Region ',num2str(k)));
            if ~isempty(neighborsavefolder)
                writetable(StatsTable,fullfile(neighborsavefolder,'/','AllNeighborStepInfo.xlsx'),'Sheet',strcat('Region ',num2str(k)));
            end
        end
    else
        skipped = skipped + 1;
    end
end
end


%% AND EXTENDED FUNCTIONS
function StatsArr = FillStatsArray(StatsArr, fullname, row, framerate, options)
i = row;
[~,ONsteps,OFFsteps,dwells,dwells_for,dwells_back] = StepandDwellHist_subfn(fullname, 0, framerate, options);

% Put whatever individual stats you want, but make sure to change it in
% FillNeighborStatsArray as well
StatsArr(i,1) = length(ONsteps(ONsteps>0)); StatsArr(i,2) = mean(ONsteps(ONsteps>0),'omitnan'); StatsArr(i,3) = std(ONsteps(ONsteps>0),'omitnan');
StatsArr(i,4) = length(ONsteps(ONsteps<0)); StatsArr(i,5) = mean(ONsteps(ONsteps<0),'omitnan'); StatsArr(i,6) = std(ONsteps(ONsteps<0),'omitnan');
StatsArr(i,7) = length(OFFsteps); StatsArr(i,8) = mean(abs(OFFsteps)); StatsArr(i,9) = std(abs(OFFsteps));
StatsArr(i,10) = length(dwells); StatsArr(i,11) = length(dwells_for); StatsArr(i,12) = mean(dwells_for);

steptrace = load(fullname);
trace = steptrace.data;
if isfield(trace,'trace') && isfield(trace,'trace_yx')
    data = trace.trace;
else
    data = [];
end
StatsArr(i,13) = size(data,1);


function NStatsArr = FillNeighborStatsArray(NStatsArr, fullname, row, framerate)
i = row;

steptrace = load(fullname);
data = steptrace.data;

% Basically do everything fNeighborlyRegions does but keep everything
% separate
if isfield(data,'trace')
    trace = data.trace;
    trace_yx = data.trace_yx;
    NearNeighborRegions = data.NeighborlyRegions;

    for k=1:length(NearNeighborRegions)
        % Steps
        [ONsteps, ~] = add_to_list_6col_steps_v3(trace,NearNeighborRegions{k},0);
        [OFFsteps, ~] = add_to_list_6col_steps_v3(trace_yx,NearNeighborRegions{k},0);

        % Dwells
        dwells = [];
        dwells_for = [];
        mat = add_to_list_6col_dwells_v3(trace,NearNeighborRegions{k},framerate,0);
        if ~isempty(mat) %check that dwells were found (JS Edit 220310)
            % All Dwells
            dwells = mat(:,3); %dwell = mat(:,3)
            % Forward and Backward dwells
            [dwells_for,dwells_back] = add_to_list_6col_dwells_for_back_v3(trace,NearNeighborRegions{k},framerate);
        end
        
        % Put whatever individual stats you want, but make sure to change it in
        % FillStatsArray as well
        NStatsArr(i,1,k) = length(ONsteps(ONsteps>0)); NStatsArr(i,2,k) = mean(ONsteps(ONsteps>0),'omitnan'); NStatsArr(i,3,k) = std(ONsteps(ONsteps>0),'omitnan');
        NStatsArr(i,4,k) = length(ONsteps(ONsteps<0)); NStatsArr(i,5,k) = mean(ONsteps(ONsteps<0),'omitnan'); NStatsArr(i,6,k) = std(ONsteps(ONsteps<0),'omitnan');
        NStatsArr(i,7,k) = length(OFFsteps); NStatsArr(i,8,k) = mean(abs(OFFsteps)); NStatsArr(i,9,k) = std(abs(OFFsteps));
        NStatsArr(i,10,k) = length(dwells); NStatsArr(i,11,k) = length(dwells_for); NStatsArr(i,12,k) = mean(dwells_for);
        NStatsArr(i,13,k) = length(NearNeighborRegions{k});
    end
end



function StatsTable = ConvertToTable(StatsArr, fnames)
% Make a table to save
StatsTable = table;
StatsTable.Name = fnames;
StatsTable.NDataPts = StatsArr(:,13);
StatsTable.ONNum = StatsArr(:,1);
StatsTable.ONMean = StatsArr(:,2);
StatsTable.ONStd = StatsArr(:,3);
StatsTable.BACKNum = StatsArr(:,4);
StatsTable.BACKMean = StatsArr(:,5);
StatsTable.BACKStd = StatsArr(:,6);
StatsTable.OFFNum = StatsArr(:,7);
StatsTable.OFFMean = StatsArr(:,8);
StatsTable.OFFStd = StatsArr(:,9);
StatsTable.AllDwells = StatsArr(:,10);
StatsTable.FORDwells = StatsArr(:,11);
StatsTable.DwellTime = StatsArr(:,12);