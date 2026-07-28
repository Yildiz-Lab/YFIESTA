function CollNeighbors = fNeighborBehavior(remove_similar, actualdir)
% Gather information about the position of neighbors relative to the trace
% and some of its dynamics
global Config

if nargin < 1
    remove_similar = 1;
end
if nargin < 2
    actualdir = uigetdir();
end
f = dir(fullfile(actualdir,'*.mat'));
fnum = length(f);

% 1 x number output args
% CollNeighbors = zeros(1,2);
CollNeighbors = nan(0,5);

for i=1:fnum
    fname = f(i).name;
    steptrace = load(strcat(actualdir,'\',fname));
    data = steptrace.data;
    if ~isfield(data,'trace') || ~isfield(data,'trace_yx')
        fnum = fnum - 1;
    else
    trace = data.trace;
    trace_yx = data.trace_yx;
    neighbors = data.neighbors;
    
    lastidx = size(CollNeighbors,1);
    
    running_MSD = cell(1,75);
    tau = 1:length(running_MSD);
    
    for n = 1:length(neighbors)
        ndata = neighbors{n}; %[rel frame, rel parallel dir, rel transverse dir]
        
        if length(ndata) > 20
        max_delta_parallel = max(ndata(:,2)) - min(ndata(:,2));
        max_delta_transverse = max(ndata(:,3)) - min(ndata(:,3));
        MSD = mean( (ndata(1:end-1,2) - ndata(2:end,2)).^2 + (ndata(1:end-1,3) - ndata(2:end,3)).^2, 'omitnan');
        
        % THIS WHOLE SECTION ONLY NECESSARY FOR MSD
            % find what preaveraging was done by looking for unique positions and
            % then take the effective number. We don't want to misdirectingly use
            % this again for MSD.
            
        % get positions of unique elements
        [~,posel] = unique(neighbors{n}(:,2:3),'rows');
        posel = sort(posel);
        preaveraged_intervals = posel(2:end) - posel(1:end-1);
        preaveraged_intervals = [preaveraged_intervals; length(neighbors{n})+1 - posel(end)];
        % chop off things that aren't at the max interval, just to make the time counting easier.
        framediff = max(preaveraged_intervals);
        posel = posel(preaveraged_intervals == framediff);
        % now do sliding window average scaling by tau
        for k = 1:length(running_MSD)
            t = tau(k);
            if length(ndata) > t*framediff
                running_MSD{t} = [running_MSD{t}; (ndata(posel(1):t*framediff:posel(end)-framediff,2) - ndata(posel(1)+framediff:t*framediff:posel(end),2)).^2 + (ndata(posel(1):t*framediff:posel(end)-framediff,3) - ndata(posel(1)+framediff:t*framediff:posel(end),3)).^2];
            end
        end
        % END OF MSD MADNESS
        
        % where is the molecule (begin, middle, end) +/- 150 nm
        relative_parallel_position = mean(ndata(:,2), 'omitnan');
        
        % to shadow whether this neighbor should be counted or not, we will
        % just translate down 200 nm (limit of "interaction") and make sure
        % the neighbor exists to the upper left of that curve
        
        % make an array that will compare the vector to connect every neighbor point
        % to the trace in the parallel direction
        ndata = ndata(~isnan(ndata(:,1)),:);
        tx = zeros(size(ndata,1),2);
        for j = 1:size(tx,1)
            td = ndata(j,1);
            if td < 1
                td = 1;
            elseif td > length(trace(:,1))-1
                td = length(trace(:,1))-1;
            end
            tx(j,1) = ndata(j,1) - td;
            tx(j,2) = ndata(j,2) - trace(td,3);
        end
        % basically check that it is in the second quadrant compared to its
        % closest values. Thus t (the x axis) must be negative and the
        % trace position x (the y axis) must be positive in at least one
        % spot
        if any( (tx(:,1) <= 0) .* (tx(:,2) > -200) )
            if relative_parallel_position < min(trace(:,3)) + 100  %if it is within or less than _ nm of the min
                tracepos = -1;
            elseif relative_parallel_position > max(trace(:,3)) - 100 %if it is within or more than _ nm of the max 
                tracepos = 1;
            else
                tracepos = 0;
            end
        else
            tracepos = NaN;
        end
        CollNeighbors = [CollNeighbors; tracepos, MSD, relative_parallel_position, max_delta_parallel, length(ndata(:,1))];
        if any(abs(CollNeighbors(lastidx+1:end-1,3) - relative_parallel_position) < 150) && remove_similar %likely similar neighbors, so remove the later ones
            CollNeighbors(end,:) = [];
        end
        end % of using this neighbor for statistics

    end
    
    end
end

% Print out some stats
size(CollNeighbors,1)

% Where are the neighbors?
fprintf(strcat("Beginning/Attachment Events: ", num2str(length(find(CollNeighbors(:,1) == -1))), "\n"))
fprintf(strcat("Ending/Detachment Events: ", num2str(length(find(CollNeighbors(:,1) == 1))), "\n"))
fprintf(strcat("Middle/Remaining Events: ", num2str(length(find(CollNeighbors(:,1) == 0))), "\n"))

% What are the MSDs of the neighbors?
fprintf(strcat("Mean Square Displacement (nm^2) / frame: ", num2str(mean(CollNeighbors(:,2)./CollNeighbors(:,5))), "\n"))

% How much do they move?
fprintf(strcat("Maximum Overall Displacement (nm) / frame: ", num2str(mean(CollNeighbors(:,4)./CollNeighbors(:,5))), "\n"))
fprintf(strcat("Maximum Overall Displacement (nm) / frame Weighted Average: ", num2str(sum(CollNeighbors(:,4)/sum(CollNeighbors(:,5)))), "\n"))

% What are their lifetimes?
framediff = framediff * Config.Time(1) / 1000;
figure()
for t = 1:length(tau)
    MSD(t) = mean(running_MSD{t},'omitnan');
    MSD_std(t) = std(running_MSD{t},'omitnan');
end
errorbar(tau*framediff, MSD, MSD_std,'o')
ylabel('nm^2')
xlabel('s')

mdl_MSD = fittype('4*D*x^a','indep','x');
mdlt1_MSD = fittype('4*D*x','indep','x');
mdl = fit(framediff*tau', MSD', mdl_MSD, 'start',[50.,1.],'weight',1./MSD_std');
mdlt1 = fit(framediff*tau', MSD', mdlt1_MSD, 'start',[50.],'weight',1./MSD_std');
hold on
plot(tau*framediff, mdl(tau*framediff), 'DisplayName',strcat("D = ", num2str(round(mdl.D,1)), ", alpha = ", num2str(round(mdl.a,2))))
plot(tau*framediff, mdlt1(tau*framediff), 'DisplayName',strcat("D = ", num2str(round(mdlt1.D,1)), ", alpha = 1"))
legend()