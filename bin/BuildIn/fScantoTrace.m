function analyzedMolecules = fScantoTrace(hMainGui)

% JS 2023/03/02
% Function to do neighbor statistics but instead on all scan region
% molecules within the region

global Config
global Stack
global Molecule

analyzedMolecules = [];

% find current channel selected
for ch = 1:length(hMainGui.ToolBar.ToolChannels)
    if strcmp(hMainGui.ToolBar.ToolChannels(ch).State, 'on')
        break
    end
end

fframe = size(Stack{1},3);

nRegion=length(hMainGui.Region);

% Store data in a separate folder for safekeeping
FFolderName = fullfile(Config.Directory{1}, strcat(Config.StackName{1}(1:end-6),'_regionscan'));
if ~isfolder(FFolderName)
    mkdir(FFolderName);
end

for j = 1:nRegion %go through each region selected
    if isfield(hMainGui.Region(j),'ScanData') % only can use scans
        fprintf(strcat("Scan to Trace on Region ", num2str(j), " \n"))
        ScanData = hMainGui.Region(j).ScanData;
        % make a FIONA-like path/trace based off this Scan Data completely
        % in the last possible frame (this way we count the statistics)
        
        % Analogous to PathData in PathStatsGui, we must make a 3xN matrix with t,X,Y
        % positions
        ScanResults = zeros(length(ScanData.X),2);
        frame1 = fframe - length(ScanData.X) + 1;
        ScanResults(:,1) = frame1:fframe;
        ScanResults(:,2) = ScanData.X';
        ScanResults(:,3) = ScanData.Y';
%         Results(:,4) = nan(length(ScanData.X),1);
        
%         [param3,~] = PathFitCurved(Results(:,2:4),4);
%         InterpPath = EvalCurvedPath(param3,InterpResults)
%         neighbors = findNeighbors(InterpPath(:,1:2));
        
        % interpolate between each data point by 1 nm
        Delta = ScanResults(2:end,2:3)-ScanResults(1:end-1,2:3);
        InterpPath = ScanResults(1,2:3); % initialize with first point
        for i = 1:size(Delta,1)
            dstep = round(Config.PixSize*sqrt(Delta(i,1)^2 + Delta(i,2)^2)); % approximately 1 nm per step
            step = Delta(i,1:2)/dstep;
            xi = ScanResults(i,2)*ones(1,dstep) + step(1)*(1:dstep);
            yi = ScanResults(i,3)*ones(1,dstep) + step(2)*(1:dstep);
            InterpPath = [InterpPath; xi' yi'];
        end
        InterpPath = Config.PixSize*InterpPath; % convert to nm at the end is fine
        
        % copied from PathStatsGui for Neighbors (maybe should make
        % separate function)
        neighbors = findNeighborsScan(Config.PixSize*[hMainGui.Region(j).X', hMainGui.Region(j).Y'], ch); %round(Config.PixSize*ScanData.ScanSize)
        neighbor_txy = cell(length(neighbors),1);
        neighbor_exist_thresh = 10;

        for m = length(neighbors):-1:1
            Res = Molecule(neighbors(m)).Results;
            if size(Res,1) < neighbor_exist_thresh
                neighbor_txy(m) = [];
            else
            txy_reorient = nan( Res(end,1)-Res(1,1) ,3);
            B = [ScanResults(end,2) - ScanResults(1,2), ScanResults(end,3) - ScanResults(1,3),0]; %for cross product
            for p = 1:size(Res,1)
    %             % Plot on direct path of molecule
    %             [npos, nidx] = min( sqrt( (PathStats(n).PathData(:,1) - Res(p,3)).^2 + (PathStats(n).PathData(:,2) - Res(p,4)).^2) );
    %             mpos = norm (PathStats(n).PathData(nidx,1:2) - PathStats(n).PathData(1,1:2));
                % Plot on interpolated path
                [npos, nidx] = min( sqrt( (InterpPath(:,1) - Res(p,3)).^2 + (InterpPath(:,2) - Res(p,4)).^2) );
                mpos = norm (InterpPath(nidx,1:2) - InterpPath(1,1:2));

                % to get correct sign off axis, use a cross product on vectors
                % for overall direction
                C = [Res(p,3) - InterpPath(nidx,1), Res(p,4) - InterpPath(nidx,2),0]; % to nearest point
                D = [Res(p,3) - ScanResults(1,2), Res(p,4) - ScanResults(1,3),0]; % to beginning
                BXC = cross(B,C); npos = npos * sign(BXC(3));
                theta = acos(dot(B,D)/norm(B)/norm(D)); mpos = mpos * sign(cos(theta));

                txy_reorient(Res(p,1)-Res(1,1)+1,:) = [Res(p,1) - frame1 + 1, mpos, npos];
            end
            neighbor_txy{m,1} = txy_reorient;
            analyzedMolecules = [analyzedMolecules, Molecule(neighbors(m)).Name];
            end
        end
        
        % Way to save a data trace, though trace data is basically meaningless
        % now. We just want it to do the statistics for the neighbors :)
        
        fname = fullfile(FFolderName,strcat('region',num2str(j)));
        xynew = zeros(length(ScanData.X),2);
        xynew(2:end,1) = Config.PixSize*sqrt(sum((ScanResults(2:end,:)-ScanResults(1:end-1,:)).^2,2));
%         if ~isfile(strcat(fname,'_fiona.mat')) % Don't overwrite existing file
            data.xy = xynew;
            data.yx = xynew(:,[2,1]);
            data.neighbors = neighbor_txy;
            % hardcode in the trace
            data.trace = [xynew, xynew, ones(size(xynew,1)), zeros(size(xynew,1))];
            data.trace_yx = [xynew(:,[2,1]), xynew(:,[2,1]), ones(size(xynew,1)), zeros(size(xynew,1))];
            save(strcat(fname,'_fiona.mat'),'data');
%         end
        
    end
    
end

function nbs = findNeighborsScan(MolXY, CheckChannel, ScanSize)
global Molecule;

if nargin<3
    ScanSize=100; % in nm, make it so that they share a point within a pixel
end

% Copied from fRightPanel NewScan
nX=MolXY(:,1);
nY=MolXY(:,2);
d=[0; cumsum(sqrt((nX(2:end)-nX(1:end-1)).^2 + (nY(2:end)-nY(1:end-1)).^2))];
dt=max(d)/round(max(d));
id=(0:round(max(d)))'*dt;
scan_length=length(id);
idx = nearestpoint(id,d);
X=zeros(scan_length,1);
Y=zeros(scan_length,1);
dis = id-d(idx);
dis(1)=0;
dis(end)=0;
X(dis==0) = nX(idx(dis==0));
Y(dis==0) = nY(idx(dis==0));
X(dis>0) = nX(idx(dis>0))+(nX(idx(dis>0)+1)-nX(idx(dis>0)))./(d(idx(dis>0)+1)-d(idx(dis>0))).*dis(dis>0);
Y(dis>0) = nY(idx(dis>0))+(nY(idx(dis>0)+1)-nY(idx(dis>0)))./(d(idx(dis>0)+1)-d(idx(dis>0))).*dis(dis>0);
X(dis<0) = nX(idx(dis<0))+(nX(idx(dis<0)-1)-nX(idx(dis<0)))./(d(idx(dis<0)-1)-d(idx(dis<0))).*dis(dis<0);
Y(dis<0) = nY(idx(dis<0))+(nY(idx(dis<0)-1)-nY(idx(dis<0)))./(d(idx(dis<0)-1)-d(idx(dis<0))).*dis(dis<0);
iX=zeros(2*ScanSize+1,scan_length);
iY=zeros(2*ScanSize+1,scan_length);
n=zeros(scan_length,3);
for i=1:length(X)
    if i==1   
        v=[X(i+1)-X(i) Y(i+1)-Y(i) 0];
        n(i,:)=[v(2) -v(1) 0]/norm(v); 
    elseif i==length(X)
        v=[X(i)-X(i-1) Y(i)-Y(i-1) 0];
        n(i,:)=[v(2) -v(1) 0]/norm(v);
    else
        v1=[X(i+1)-X(i) Y(i+1)-Y(i) 0];
        v2=-[X(i)-X(i-1) Y(i)-Y(i-1) 0];
        n(i,:)=v1/norm(v1)+v2/norm(v2); 
        if norm(n(i,:))==0
            n(i,:)=[v1(2) -v1(1) 0]/norm(v1);
        else
            n(i,:)=n(i,:)/norm(n(i,:));
        end
        z=cross(v1,n(i,:));
        if z(3)>0
            n(i,:)=-n(i,:);
        end
    end
    iX(:,i)=linspace(X(i)+ScanSize*n(i,1),X(i)-ScanSize*n(i,1),2*ScanSize+1)';
    iY(:,i)=linspace(Y(i)+ScanSize*n(i,2),Y(i)-ScanSize*n(i,2),2*ScanSize+1)';
end

iX = reshape(iX,numel(iX),1);
iY = reshape(iY,numel(iY),1);
% find boundary for our polygon
bd = boundary(iX,iY);

% check if in polygon for all molecules in other channel (rudimentary, but
% fair)
nbs = [];
for j = 1:length(Molecule)
    if Molecule(j).Channel == CheckChannel %only do if in opposite channel
        inbd = inpolygon(Molecule(j).Results(:,3), Molecule(j).Results(:,4), iX(bd), iY(bd));
        if any(inbd > 0)
            nbs = [nbs, j];
        end
    end
end


function CollNeighbors = fNeighborBehaviorScan(actualdir)
% Gather information about the position of neighbors relative to the trace
% and some of its dynamics

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
    
    for n = 1:length(neighbors)
        ndata = neighbors{n}; %[rel frame, rel parallel dir, rel transverse dir]
        max_delta_parallel = max(ndata(:,2)) - min(ndata(:,2));
        max_delta_transverse = max(ndata(:,3)) - min(ndata(:,3));
        MSD = mean( (ndata(1:end-1,2) - ndata(2:end,2)).^2 + (ndata(1:end-1,3) - ndata(2:end,3)).^2, 'omitnan');
        
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
        CollNeighbors = [CollNeighbors; tracepos, MSD, relative_parallel_position, max_delta_parallel, length(ndata(:,1))];
%         if any(abs(CollNeighbors(lastidx+1:end-1,3) - relative_parallel_position) < 150) %likely similar neighbors, so remove the later ones
%             CollNeighbors(end,:) = [];
%         end
    end
    
    end
end

% Print out some stats

% What are the MSDs of the neighbors?
fprintf(strcat("Mean Square Displacement (nm^2) / frame: ", num2str(mean(CollNeighbors(:,2)./CollNeighbors(:,5))), "\n"))

% How much do they move?
fprintf(strcat("Maximum Overall Displacement (nm) / frame: ", num2str(mean(CollNeighbors(:,4)./CollNeighbors(:,5))), "\n"))
CollNeighbors(:,4)./CollNeighbors(:,5)


% What are their lifetimes?


