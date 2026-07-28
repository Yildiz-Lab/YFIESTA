function fGenerateFIONA4Neighbors(minj, nbhj, PathStats_n, txy_reorient_original, InterpPath, fname, options)
% minj - original molecule number, index version
% nbhj - neighbor molecule number, index version
% txy_reorient_original - Nx5 array holding time, distance in x, and distance in y
% InterpPath - Interpolated path according to original molecule
% fname - root filename
% options - passed in options from fNeighborStepOptions still apply

global Molecule;

[dir, fname, ~] = fileparts(fname);

frames = Molecule(nbhj).Results(:,1); frame1 = frames(1);
xynew = txy_reorient_original(:,2:3);
idx = find(~isnan(txy_reorient_original(:,4)));
nghb_molecule_num = txy_reorient_original(idx(1),4);
time = txy_reorient_original(:,5);

% JS Edit 2024/12/01 copied straight from fPathStatsGUI > FIONA
if Molecule(nbhj).Channel > 1
    neighbors = findNeighbors(InterpPath(:,1:2), options, 1);
else
    neighbors = findNeighbors(InterpPath(:,1:2), options);
end
    
neighbor_txy = cell(length(neighbors),1);

for m = length(neighbors):-1:1
    Res = Molecule(neighbors(m)).Results;
    neighbor_footprint = [Res(1,2) - options.eExcludeTime, Res(end,2) + options.eExcludeTime];
    % if size(Res,1) < options.eExcludeTime
    %     neighbor_txy(m) = [];
    if Res(end,2)-Res(1,2) < options.ExistThresh % time is too short
        neighbor_txy(m) = []; %not a long enough lifetime
     %   Exclude Time       &&         (  completely before || completely after  )
    elseif options.ExcludeTime && ( all( neighbor_footprint < Molecule(minj).Results(1,2)) ||  all( neighbor_footprint > Molecule(minj).Results(end,2)) )
        neighbor_txy(m) = []; %our neighbor is too far before the molecule or too far after
    else
    % txy_reorient = nan( Res(end,1)-Res(1,1) , 3);
    txy_reorient = nan( Res(end,1)-Res(1,1) , 5);
    B = [PathStats_n.Results(end,3) - PathStats_n.Results(1,3), PathStats_n.Results(end,4) - PathStats_n.Results(1,4),0]; %for cross product
    for p = 1:size(Res,1)
%             % Plot on direct path of molecule
%             [npos, nidx] = min( sqrt( (PathStats_n.PathData(:,1) - Res(p,3)).^2 + (PathStats_n.PathData(:,2) - Res(p,4)).^2) );
%             mpos = norm (PathStats_n.PathData(nidx,1:2) - PathStats_n.PathData(1,1:2));
        % Plot on interpolated path
        [npos, nidx] = min( sqrt( (InterpPath(:,1) - Res(p,3)).^2 + (InterpPath(:,2) - Res(p,4)).^2) );
        mpos = norm (InterpPath(nidx,1:2) - PathStats_n.PathData(1,1:2));
        
        % to get correct sign off axis, use a cross product on vectors
        % for overall direction
        C = [Res(p,3) - InterpPath(nidx,1), Res(p,4) - InterpPath(nidx,2),0]; % to nearest point
        D = [Res(p,3) - PathStats_n.Results(1,3), Res(p,4) - PathStats_n.Results(1,4),0]; % to beginning
        BXC = cross(B,C); npos = npos * sign(BXC(3));
        theta = acos(dot(B,D)/norm(B)/norm(D)); mpos = mpos * sign(cos(theta));
        
        % txy_reorient(Res(p,1)-Res(1,1)+1,:) = [Res(p,1) - frame1 + 1, mpos, npos];
        % JS Edit 2024/11/27 to get the molecule and molecule time data
        % saved
        txy_reorient(Res(p,1)-Res(1,1)+1,:) = [Res(p,1) - frame1 + 1, mpos, npos, str2double(Molecule(neighbors(m)).Name(10:end)), Res(p,2)];
    end
    neighbor_txy{m,1} = txy_reorient;

    end
    
end

fname_new = fullfile(dir, strcat(fname, '_nbh_', num2str(nghb_molecule_num)));

% New way to save
if ~isfile(strcat(fname_new,'_fiona.mat')) % Don't overwrite existing file
data.xy = xynew;
data.yx = xynew(:,[2,1]);
data.neighbors = neighbor_txy;
data.time = time;
save(strcat(fname_new,'_fiona.mat'),'data');
end


%JS Edit 2022/09/26 for finding nearby molecules from fPathStatsGUI
function nbs = findNeighbors(MolXY, options, CheckChannel)
global Molecule;

if nargin < 3
    CheckChannel = 2;
end

% Copied from fRightPanel NewScan
nX=MolXY(:,1);
nY=MolXY(:,2);
% ScanSize=100; % in nm, make it so that they share a point within a pixel
% ScanSize=10; % in nm, MINFLUX we can make this a bit closer
ScanSize=max(options.XA(2), options.XB(2)); %JS Edit 2024/11/27
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