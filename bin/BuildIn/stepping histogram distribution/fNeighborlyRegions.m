function fNeighborlyRegions(opts, framerate, actualdir, noplot, savedir)
%opts should have field mode and then any other parameters that specific
%mode might require
% Opts should include XA, XB, YA, YB

if nargin < 3
    actualdir = uigetdir();
end
f = dir(fullfile(actualdir,'*.mat'));
fnum = length(f);
%% Choose your regions in x and y to consider ranges of close interactions.
% Regions are arrays defining distance away from mean MAP location.
% Example: [0,50,100] means Region1 [0,50), Region2 [50,100), Region3 [100,Inf)
% b = before meaning as the motor approaches the MAP mean position
% a = after meaning after the motor has passed the MAP mean position
% after is to exploit possible asymmetry in the MAPs behavior. If you want
% symmetry, just call them equivalent

% deal with all the cases of no parameters passed, here are defaults
% assume symmetry for any cases less than all

if nargin < 4
    noplot = 0;
end
if nargin < 5
    savedir = [];
end
xb = opts.XB; xa = opts.XA; yb = opts.YB; ya = opts.YA;

% Storage for statistics in each region
do = struct();
switch(opts.mode)
    case 'steps'
        do.steps = 1;
        RegionStepStats = cell(1,length(xb));
        RegionOffStepStats = cell(1,length(xb));
        RegionDwellStats = cell(1,length(xb));
        RegionDwellForStats = cell(1,length(xb));
        RegionDwellBackStats = cell(1,length(xb));
    case 'pause'
        do.pause = 1;
        RegionPauseFreqStats = cell(1,length(xb));
        RegionPauseFracStats = cell(1,length(xb));
        RegionPauseStats = cell(1,length(xb));
        RegionPauseSideStats = cell(1,length(xb));
        RegionPauseBackStats = cell(1,length(xb));
    case 'regions'
        fprintf("Just finding regions")
end

% A simple count to keep track of how many neighbors we have
ncount=0;

%% Now load each MATLAB file and determine neighbor regions
for i=1:fnum
    fname = f(i).name;
    steptrace = load(strcat(actualdir,'\',fname));
    data = steptrace.data;
    if ~isfield(data,'trace') || ~isfield(data,'trace_yx')
        fnum = fnum - 1;
    else
    trace = data.trace;
    trace_yx = data.trace_yx;

    RM = zeros(length(data.neighbors),6);
    Regions = cell(1,length(xb));
    % RM [time start, time end, meanx, stdx, meany, stdy]
    ncount = ncount + length(data.neighbors);
    for n=1:length(data.neighbors)
%         rsummod = data.neighbors{n}(1:sumidx:end,:);
%         RM(n,1) = rsummod(1,1); RM(n,2) = rsummod(end,1);
%         RM(n,3) = mean(rsummod(:,2)); RM(n,4) = std(rsummod(:,2));
%         RM(n,5) = mean(rsummod(:,3)); RM(n,6) = std(rsummod(:,3));
        r = data.neighbors{n};
        RM(n,1) = r(1,1); RM(n,2) = r(end,1);
        RM(n,3) = mean(r(:,2),'omitnan'); RM(n,4) = std(r(:,2),'omitnan');
        RM(n,5) = mean(r(:,3),'omitnan'); RM(n,6) = std(r(:,3),'omitnan');
    end

    % Now RM has all the information to go through, but most importantly
    % t,tend,meanx,stdx,meany,stdy for each region defined (think like box border)     
    % it is 3D trace: t, on-axis(x), off-axis(y)
    % and is it contained in box [t,tend; xlimits; ylimits]
    for b=2:length(xb)
        for n=1:length(data.neighbors)
            tb = RM(n,1); xc = RM(n,3); yc = RM(n,5);
            k = inpolygon(trace(:,1), trace(:,2), xc + [-xb(b),-xb(b),xa(b),xa(b),-xb(b)], yc + [-yb(b),ya(b),-yb(b),ya(b),-yb(b)]);
            [t,~] = find(k == 1);
            %these time (actually index) points should now be considered in this region
            % and only if they are greater than the time when the neighbor
            % first appeared
            t = t(t>tb);
            % In the end, union the sets
            Regions{b-1} = unique([Regions{b-1} t']);
        end

        % those that haven't been matched by the end should be considered
        % in the final region
        if b == length(xb)
            notnan = 1:length(trace);
            Regions{b} = notnan(~isnan(trace(:,1)'));
        end

    end

    % Finally, if we have a previous region, and we want to exclude
    % this behavior (like rings rather than a convex circle), then you
    % should make sure these t points don't exist in the previous
    % region. This can be optioned out very easily by another passed
    % param
    for m=1:length(xb)
        for n=m-1:-1:1
            RB = Regions{m};
            Regions{m} = RB(~ismember(RB,Regions{n}));
        end
    end

    % we now have time points that are in user defined region. Now, there
    % remain two issues possible:
        %  - if something is on the borderline, then it will alternate between
        %  two regions. We can fix this by picking the step that is most
        %  associated with this region
        %  - how do we pass this data in the future? Let's keep it as time
        %  regions which will allow us to index later

    % Postprocessing Check with steps in dominant axis (x)
    % get changepoints for step transitions and include NaNs again
    NearNeighborRegions = cell(1,length(xb));
    [idx,~] = find(trace(:,5)==1);
    idx = [1,idx',length(trace)];
    for j = 1:length(idx)-1
        nummember = -1;
        %for b = length(xb):-1:1
        for b = 1:length(xb)
            tofind = idx(j):idx(j+1)-1;
            b;
            ismember(tofind,Regions{b});
            num = length(tofind(ismember(tofind,Regions{b})==1));
            if num > nummember
                nummember = num;
                regionplace = b;
            end
        end
        NearNeighborRegions{regionplace} = [NearNeighborRegions{regionplace} tofind];
    end

    % Now these regions are prepared to ship off to data processing. Add
    % this to the files for storage. But meanwhile we can begin compiling
    % information for statistics.
    data.NeighborlyRegions = NearNeighborRegions;
    save(fullfile(actualdir,fname),'data');

    %% Now starts storage options

    if isfield(do,'steps')
        for k=1:length(NearNeighborRegions)
            % Steps
            [on_steps, off_steps] = add_to_list_6col_steps_v3(trace,NearNeighborRegions{k},0);
            RegionStepStats{k} = [RegionStepStats{k}; on_steps'];

            [on_steps, off_steps] = add_to_list_6col_steps_v3(trace_yx,NearNeighborRegions{k},0);
            RegionOffStepStats{k} = [RegionOffStepStats{k}; on_steps'];

            % Dwells
            mat = add_to_list_6col_dwells_v3(trace,NearNeighborRegions{k},framerate,0);
            if ~isempty(mat) %check that dwells were found (JS Edit 220310)
                % All Dwells
                RegionDwellStats{k} = [RegionDwellStats{k}; mat(:,3)]; %dwell = mat(:,3)

                % Forward and Backward dwells
                [forward,backward] = add_to_list_6col_dwells_for_back_v3(trace,NearNeighborRegions{k},framerate);
                RegionDwellForStats{k} = [RegionDwellForStats{k}; forward];
                RegionDwellBackStats{k} = [RegionDwellBackStats{k}; backward];
            end
        end
    end

    if isfield(do,'pause')
        for k=1:length(NearNeighborRegions)
            % JS Edit 2023/01/09 Pause Stats based on DeWitt et al (2015)
            if ~isempty(NearNeighborRegions{k})
                % Note we are always truncating off the last point
                [pause_frequency, pause_fraction, ~, numpause, numside, numback] = fPauseAnalysis(trace, trace_yx, NearNeighborRegions{k}, opts.PauseThresh);
                RegionPauseFreqStats{k} = [RegionPauseFreqStats{k}; pause_frequency];
                RegionPauseFracStats{k} = [RegionPauseFracStats{k}; pause_fraction];
                RegionPauseStats{k} = [RegionPauseStats{k}; numpause];
                RegionPauseSideStats{k} = [RegionPauseSideStats{k}; numside];
                RegionPauseBackStats{k} = [RegionPauseBackStats{k}; numback];
            end
        end
    end

    end %of isfield data
end

%% Now compile differently

if isfield(do,'steps')

    if ~isempty(RegionStepStats{1}) && ~noplot
        for k=1:length(xb)
            if ~isempty(savedir)
                PlotStepStats(fnum,RegionStepStats{k},RegionOffStepStats{k},RegionDwellStats{k},RegionDwellForStats{k},RegionDwellBackStats{k},fullfile(savedir,strcat('Region',num2str(k))))
            else
                PlotStepStats(fnum,RegionStepStats{k},RegionOffStepStats{k},RegionDwellStats{k},RegionDwellForStats{k},RegionDwellBackStats{k})
            end
        end
    end
    fprintf(strcat("Number of Neighbors Used in Stats: ", num2str(ncount), "\n"))

    %% Make a summary of all the regions information

    % only do summary plot if there are neighbors to do so.
    % At the moment, if one Region is empty, we throw it away.
    if ~noplot && ~isempty(RegionStepStats{1})
        h =  findobj('type','figure');
        n = length(h);
        figure(n+1)
        % Back Steps
        subplot(2,2,4)
        for k=1:length(xb)
            Nsteps(k) = length(RegionStepStats{k});
            forsteps(k) = length(RegionStepStats{k}(RegionStepStats{k} > 0));
            backsteps(k) = length(RegionStepStats{k}(RegionStepStats{k} < 0));
            sidesteps(k) = length(RegionOffStepStats{k});
            obj0 = cdfplot(RegionDwellForStats{k});
            mdl_gamma_cdf = fittype('real(gammainc(k*x,2))','indep','x');
            X = obj0.XData(2:end-1); % get rid of infs
            Y = obj0.YData(2:end-1);
            fittedmdl = fit(X',Y',mdl_gamma_cdf,'start',[3.]);
            cint(1:2,k) = confint(fittedmdl);
            krate(k) = fittedmdl.k;
        end

        subplot(2,2,1)
        b1 = bar(1:length(backsteps),backsteps./Nsteps);
        xtips1 = b1.XEndPoints;
        ytips1 = b1.YEndPoints;
        labels1 = cell(1,length(xb));
        for i = 1:length(xb)
            labels1{i} = strcat(num2str(backsteps(i)),'/',num2str(Nsteps(i)));
        end
        text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
        title("Percentage Back Steps")

        subplot(2,2,2)
        b2 = bar(1:length(sidesteps),sidesteps./Nsteps);
        xtips2 = b2.XEndPoints;
        ytips2 = b2.YEndPoints;
        labels2 = cell(1,length(xb));
        for i = 1:length(xb)
            labels2{i} = strcat(num2str(sidesteps(i)),'/',num2str(Nsteps(i)));
        end
        text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
        title("Percentage Side Steps")

        subplot(2,2,3)
        b3 = bar(krate);
        hold on
        cerr = (cint(2,:)-cint(1,:))/2;
        errorbar(b3.XEndPoints,b3.YEndPoints,cerr,'k.','CapSize',22,'LineWidth',2)
        xtips3 = b3.XEndPoints;
        ytips3 = b3.YEndPoints-cerr-0.2; %extra 0.2 for font size
        labels3 = string(krate);
        text(xtips3,ytips3,labels3,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
        xlabel("Region")
        title("Stepping rate")

    end
end

if isfield(do,'pause')
    for k=1:length(NearNeighborRegions)
        fprintf(strcat("Region ",num2str(k)," pause number: ", num2str(sum(RegionPauseStats{k})), "\n"))
        fprintf(strcat("Region ",num2str(k)," pause density: ", num2str(mean(RegionPauseFreqStats{k})*1000), "\n"))
        fprintf(strcat("Region ",num2str(k)," pause fraction: ", num2str(mean(RegionPauseFracStats{k})), "\n"))
        fprintf(strcat("Region ",num2str(k)," pause-side fraction: ", num2str(sum(RegionPauseSideStats{k})/sum(RegionPauseStats{k})), " (", num2str(sum(RegionPauseSideStats{k})), "/", num2str(sum(RegionPauseStats{k})), ")", "\n"))
        fprintf(strcat("Region ",num2str(k)," pause-back fraction: ", num2str(sum(RegionPauseBackStats{k})/sum(RegionPauseStats{k})), " (", num2str(sum(RegionPauseBackStats{k})), "/", num2str(sum(RegionPauseStats{k})), ")", "\n"))
    end
end

end

