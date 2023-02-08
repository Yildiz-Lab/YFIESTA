function CollNeighbors = fNeighborBehavior(actualdir)
% Gather information about the position of neighbors relative to the trace
% and some of its dynamics

if nargin < 2
    actualdir = uigetdir();
end
f = dir(fullfile(actualdir,'*.mat'));
fnum = length(f);

% 1 x number output args
% CollNeighbors = zeros(1,2);
CollNeighbors = nan(0,4);

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
    
    for n = 1:length(neighbors)
        ndata = neighbors{n}; %[rel frame, rel parallel dir, rel transverse dir]
        max_delta_parallel = max(ndata(:,2)) - min(ndata(:,2));
        max_delta_transverse = max(ndata(:,3)) - min(ndata(:,3));
        MSD = mean( (ndata(1:end-1,2) - ndata(2:end,2)).^2 + (ndata(1:end-1,3) - ndata(2:end,3)).^2, 'omitnan');
        
        % where is the molecule (begin, middle, end) +/- 150 nm
        relative_parallel_position = mean(ndata(:,2), 'omitnan');
        
        % to shadow whether this neighbor should be counted or not, we will
        % just translate down 200 nm (limit of "interaction") and make sure
        % the MAP exists to the upper left of that curve
        
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
            if relative_parallel_position < min(trace(:,3)) + 100
                tracepos = -1;
            elseif relative_parallel_position > max(trace(:,3)) - 100
                tracepos = 1;
            else
                tracepos = 0;
            end
        else
            tracepos = NaN;
        end
        CollNeighbors = [CollNeighbors; tracepos, MSD, max_delta_parallel, max(ndata(:,1)) - min(ndata(:,1))];
        
    end
    
    end
end

end

