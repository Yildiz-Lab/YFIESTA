function [Region,Alignment] = fAlignRegionsToProjection(Region,Projection,MaxShift)
%FALIGNREGIONSTOPROJECTION Register scan centerlines by one x/y translation.
% Coordinates and displacements are pixels; +X is right, +Y is down.
% Match equally weighted centerlines to multiscale bright image ridges,
% suppressing round spots. Search all integer shifts, then refine subpixel.
% Region shapes, scan widths, and measured lengths are preserved. Area masks
% are regenerated for the target image. No image or source file is modified.
if nargin<3, MaxShift=min(32,floor(min(size(Projection))/4)); end
validateattributes(Projection,{'numeric'},{'real','2d','nonempty','finite'});
validateattributes(MaxShift,{'numeric'},{'real','scalar','integer','positive','finite'});
if isempty(Region)
    error('FIESTA:NoRegions','There are no regions to align.');
end
image=log1p(max(double(Projection),0));
if max(image(:))-min(image(:))<1e-10
    error('FIESTA:NoAlignmentSignal','The z-projection has no contrast for region alignment.');
end
points=[];
weights=[];
for k=1:numel(Region)
    if ~isfield(Region(k),'ScanData') || ~isstruct(Region(k).ScanData) ...
            || ~all(isfield(Region(k).ScanData,{'X','Y'}))
        error('FIESTA:NoRegionCenterline', ...
            'Alignment requires ScanData.X/Y centerlines for every loaded region (region %d).',k);
    end
    X=double(Region(k).ScanData.X(:));
    Y=double(Region(k).ScanData.Y(:));
    if numel(X)~=numel(Y) || numel(X)<2 || any(~isfinite([X;Y]))
        error('FIESTA:NoRegionCenterline','Region %d has an invalid centerline.',k);
    end
    keep=[true; hypot(diff(X),diff(Y))>0];
    X=X(keep); Y=Y(keep);
    if numel(X)<2
        error('FIESTA:NoRegionCenterline','Region %d has a zero-length centerline.',k);
    end
    d=[0; cumsum(hypot(diff(X),diff(Y)))];
    t=linspace(0,d(end),max(2,round(d(end))+1))';
    sampled=[interp1(d,X,t) interp1(d,Y,t)];
    points=[points;sampled]; %#ok<AGROW>
    weights=[weights;ones(numel(t),1)/numel(t)]; %#ok<AGROW>
end
weights=weights/sum(weights);
ridge=RidgeImage(image);
if max(ridge(:))<=0
    error('FIESTA:NoAlignmentSignal','No filament-like signal was found in the z-projection.');
end
shifts=-MaxShift:MaxShift;
scores=zeros(numel(shifts),numel(shifts));
% Keep the denominator fixed: moving samples out of view cannot improve
% a candidate merely by dropping its poorly matched points.
for row=1:numel(shifts)
    queryX=bsxfun(@plus,points(:,1),shifts);
    queryY=repmat(points(:,2)+shifts(row),1,numel(shifts));
    values=interp2(ridge,queryX,queryY,'linear',0);
    scores(row,:)=weights'*values;
end
[bestScore,best]=max(scores(:));
[row,column]=ind2sub(size(scores),best);
seed=[shifts(column) shifts(row)];
if max(abs(seed))==MaxShift
    error('FIESTA:AlignmentSearchLimit', ...
        'The best shift reaches the %g-pixel search limit; regions were not aligned.',MaxShift);
end
significance=(bestScore-mean(scores(:)))/max(std(scores(:)),eps);
if significance<3
    error('FIESTA:AmbiguousAlignment', ...
        'The z-projection does not give a distinct alignment peak; regions were not aligned.');
end
options=optimset('Display','off','TolX',1e-3,'TolFun',1e-8,'MaxIter',200);
shift=fminsearch(@(s) TranslationCost(s,ridge,points,weights,MaxShift),seed,options);
refinedScore=-TranslationCost(shift,ridge,points,weights,MaxShift);
if refinedScore<bestScore
    shift=seed;
    refinedScore=bestScore;
end
Alignment=struct('ShiftXY',shift,'ScoreBefore', ...
    -TranslationCost([0 0],ridge,points,weights,MaxShift), ...
    'ScoreAfter',refinedScore,'PeakSignificance',significance, ...
    'MaxShift',MaxShift,'ImageSize',size(Projection));
for k=1:numel(Region)
    Region(k).X=Region(k).X+shift(1);
    Region(k).Y=Region(k).Y+shift(2);
    Region(k).ScanData.X=Region(k).ScanData.X+shift(1);
    Region(k).ScanData.Y=Region(k).ScanData.Y+shift(2);
    Region(k).Area=roipoly(size(Projection,1),size(Projection,2),Region(k).X,Region(k).Y);
end

function cost=TranslationCost(shift,ridge,points,weights,MaxShift)
if any(abs(shift)>MaxShift)
    cost=1+sum(max(abs(shift)-MaxShift,0).^2);
else
    values=interp2(ridge,points(:,1)+shift(1),points(:,2)+shift(2),'linear',0);
    cost=-weights'*values;
end

function ridge=RidgeImage(image)
ridge=zeros(size(image));
for sigma=[1 2]
    x=-ceil(4*sigma):ceil(4*sigma);
    g=exp(-x.^2/(2*sigma^2)); g=g/sum(g);
    dg=-x.*g/sigma^2;
    ddg=(x.^2/sigma^4-1/sigma^2).*g;
    ddg=ddg-sum(ddg)*g; % exactly zero response to constant background
    xx=imfilter(imfilter(image,g','replicate'),ddg,'replicate')*sigma^2;
    yy=imfilter(imfilter(image,ddg','replicate'),g,'replicate')*sigma^2;
    xy=imfilter(imfilter(image,dg','replicate'),dg,'replicate')*sigma^2;
    delta=sqrt((xx-yy).^2+4*xy.^2);
    a=(xx+yy+delta)/2; b=(xx+yy-delta)/2;
    swap=abs(a)>abs(b);
    small=a; large=b;
    small(swap)=b(swap); large(swap)=a(swap);
    response=max(-large,0).*exp(-0.5*(small./(0.5*max(abs(large),eps))).^2);
    ridge=max(ridge,response);
end
scale=prctile(ridge(:),99);
if scale>0
    ridge=imgaussfilt(min(ridge/scale,1),1,'Padding','replicate');
end
