function [Scan,X,Y] = fScanGeometry(X,Y,ScanSize)
%FSCANGEOMETRY Sample a pixel-space centerline and its scan boundaries.
% Matches fRightPanel/NewScan geometry for ROI import. ScanSize is the
% half-width in pixels; the intensity scan has 2*ScanSize+1 samples.
validateattributes(X,{'numeric'},{'real','vector','finite','nonempty'});
validateattributes(Y,{'numeric'},{'real','vector','finite','nonempty'});
validateattributes(ScanSize,{'numeric'},{'real','scalar','integer','positive','finite'});
if numel(X)~=numel(Y)
    error('FIESTA:InvalidScan','Scan X and Y must have the same length.');
end
X=double(X(:));
Y=double(Y(:));
% Repeated vertices have no direction and would produce NaN normals.
keep=[true; diff(X)~=0 | diff(Y)~=0];
X=X(keep);
Y=Y(keep);
if numel(X)<2
    error('FIESTA:InvalidScan','A scan needs at least two distinct points.');
end
Scan.X=X';
Scan.Y=Y';
nX=X(:);
nY=Y(:);
d=[0; cumsum(sqrt((nX(2:end)-nX(1:end-1)).^2 + (nY(2:end)-nY(1:end-1)).^2))];
nIntervals=max(1,round(max(d)));
dt=max(d)/nIntervals;
id=(0:nIntervals)'*dt;
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
d = [0; cumsum(sqrt((X(2:end)-X(1:end-1)).^2 + (Y(2:end)-Y(1:end-1)).^2))];
lx = iX(1,:);
ux = iX(end,:);
ly = iY(1,:);
uy = iY(end,:);
Scan.InterpX=iX;
Scan.InterpY=iY;
Scan.lx=lx;
Scan.ux=ux;
Scan.ly=ly;
Scan.uy=uy;
Scan.InterpD=d;
