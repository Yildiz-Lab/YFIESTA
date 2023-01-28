function fPauseAnalysis_wrapper(directory)
%Do Pause Analysis for trajectories evaluated by path statistics
global Config;

if length(Config.Time) > 1
    framerate = Config.Time(1)/1000;
else
    framerate = Config.Time/1000;
end

if isfield(Config,'PauseThreshold') %presaved options load into GUI
    prevopts.PauseThreshold = Config.PauseThreshold;
    opts = fPauseAnalysisOpts(prevopts);
else
    opts = fPauseAnalysisOpts();
end

if isempty(opts)
    return
end

if nargin < 1
    checkdir = Config.Directory{1};
    directory = uigetdir(checkdir);
end

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

HistVals = [];
pause_frequency = [];

% I should always have an option whether you need data or not

if opts.UseNeighborRegions
    
    if isfield(Config,'NeighborRegionsOpts') %presaved options load into GUI
        prevopts = Config.NeighborRegionsOpts;
        prevopts.UseNeighborRegions = 1;
        opts2 = fNeighborRegionOptions(prevopts);
    else
        opts2 = fNeighborRegionOptions();
    end
    %Options to put in regions here, maybe if a GUI comes
%     xb = [0,200]; yb = [0,125];
%     %xb = [0,100,200]; yb = [0,200,200];
%     xa = xb; ya = yb;
    %Forward/Backward Scheme
    % xb = 1.5*[0,50,50,100,100]; yb = [0,35,35,70,70];
    % xa = 1.5*[0,0,50,50,100]; ya = yb;

    % Generate an automatic foldername that carries Neighbor Info
    totarr = [opts2.XB,opts2.YB,opts2.XA,opts2.YA];
    foldername = '[';
    for j = 1:length(totarr)
        foldername = strcat(foldername, num2str(totarr(j)), ',');
        if mod(j,length(opts2.XB)) == 0 && j ~= length(totarr)
            foldername = strcat(foldername(1:end-1), '],[');
        end
    end
    foldername = strcat(foldername(1:end-1), ']');
    foldername = fullfile(directory,foldername);
    if ~isfolder(foldername)
        mkdir(foldername)
    end

    opts2.mode = 'pause';
    opts2.PauseThresh = opts.PauseThresh;
    opts2.UseNeighborRegions = 1;
    fNeighborlyRegions(opts2,framerate,directory,0,foldername)
    
else
    for i=1:fnum
    
        if isfile(directory)
            steptrace = load(fname);
        else
            fname = f(i).name;
            steptrace = load(strcat(cd,'\',fname));
        end

        trace = steptrace.data;
        if ~isfield(trace,'trace') % I should still find a way to look at this even if there is no trace
            fnum = fnum - 1;
        else
        data = trace.trace;
        data_yx = trace.trace_yx;

        % for pause events
        % cnt_pause_events = 47; %make sure to also change in Neighbor Regions
        % Note we are always truncating off the last point in Neighbor Regions
        % since looking for transitions (which require two data points).
        % Therefore we will also truncate exactly one data point per molecule
        % to be consistent
        [pf, Values] = fPauseAnalysis(data(1:end-1,:), [], opts.PauseThresh);
        HistVals = [HistVals, Values];
        pause_frequency = [pause_frequency, pf];
        end
    end 
    
end


% To Set Pause Threshold

if opts.PauseThresh == 0
    figure()
    hh = histogram(HistVals,'BinWidth',1);
    hold on
    CSum = zeros(1,hh.NumBins);
    CSum(1) = hh.Values(1);
    for b = 2:hh.NumBins
        CSum(b) = CSum(b-1) + hh.Values(b);
    end
    Cnorm = CSum/CSum(end);
    % find tau by finding where population is at 1/e
    [~,tau] = min(abs(Cnorm - exp(-1)));
    % find actual 95% confidence interval by finding where population is at 95%
    [~,conf2] = min(abs(Cnorm - 0.95))
    plot((conf2+1)*ones(1,2),[0,max(hh.Values)],'r--');
    Config.PauseThreshold = conf2;
elseif ~opts.UseNeighborRegions
    figure()
    hh = histogram(HistVals,'BinWidth',1);
    hold on
    plot((opts.PauseThresh+1)*ones(1,2),[0,max(hh.Values)],'r--');
    fprintf(strcat("Pause frequency: ", num2str(mean(pause_frequency)), "\n"))
end

end
