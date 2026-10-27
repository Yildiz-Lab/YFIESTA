function varargout=fTraceMovieBackground(action,varargin)
%FTRACEMOVIEBACKGROUND Shared stack sampling for preview, projection and export.
switch lower(action)
    case 'validate',varargout{1}=Validate(varargin{1});
    case 'prepare',varargout{1}=Prepare(varargin{:});
    case 'frame',[varargout{1:nargout}]=Frame(varargin{:});
    case 'projection',varargout{1}=Projection(varargin{:});
end

function b=Validate(b)
if ~isstruct(b)||~isscalar(b)||~all(isfield(b,{'Stack','TimesS','PixelSizeNm','Limits'}))
    error('FIESTA:TraceMovieBackground','Background requires Stack, TimesS, PixelSizeNm and Limits.');
end
validateattributes(b.Stack,{'numeric'},{'real','nonempty','finite'});
if ndims(b.Stack)>3,error('FIESTA:TraceMovieBackground','Use one grayscale stack channel.');end
validateattributes(b.TimesS,{'numeric'},{'real','finite','vector','numel',size(b.Stack,3)});
b.TimesS=double(b.TimesS(:));
if any(diff(b.TimesS)<=0),error('FIESTA:TraceMovieBackground','Stack times must increase, in seconds.');end
validateattributes(b.PixelSizeNm,{'numeric'},{'real','finite','positive','scalar'});
validateattributes(b.Limits,{'numeric'},{'real','finite','vector','numel',2});
b.Limits=double(reshape(b.Limits,1,2));
if b.Limits(2)<=b.Limits(1),error('FIESTA:TraceMovieBackground','Background display maximum must exceed its minimum.');end
if ~isfield(b,'FrameDurationS')||isempty(b.FrameDurationS)
    b.FrameDurationS=1;if numel(b.TimesS)>1,b.FrameDurationS=median(diff(b.TimesS));end
end
validateattributes(b.FrameDurationS,{'numeric'},{'real','finite','positive','scalar'});
if ~isfield(b,'Drift'),b.Drift=[];end
if ~isempty(b.Drift)
    validateattributes(b.Drift,{'numeric'},{'real','finite'});
    if size(b.Drift,1)~=3||size(b.Drift,2)~=3||~ismember(size(b.Drift,3),[1 size(b.Stack,3)])
        error('FIESTA:TraceMovieBackground','Drift must be one 3x3 matrix or one per stack frame.');
    end
    for k=1:size(b.Drift,3)
        T=b.Drift(:,:,k);
        if abs(det(T))<1e-12||norm(T(:,3)-[0;0;1])>1e-8
            error('FIESTA:TraceMovieBackground','Stack drift transforms must be invertible affine matrices.');
        end
    end
end

function prepared=Prepare(b,plan,sizePx)
b=Validate(b);span=diff(plan.BoundsNm,1,1);
[x,y]=meshgrid(plan.BoundsNm(1,1)+((1:sizePx(2))-.5)*span(1)/sizePx(2), ...
    plan.BoundsNm(1,2)+((1:sizePx(1))-.5)*span(2)/sizePx(1));
R=plan.RotationMatrix;c=plan.CenterNm;
qx=(x*R(1,1)+y*R(1,2)+c(1))/b.PixelSizeNm;
qy=(x*R(2,1)+y*R(2,2)+c(2))/b.PixelSizeNm;
prepared=struct('Source',b,'X',qx,'Y',qy,'Size',sizePx);

function [gray,index]=Frame(prepared,time)
b=prepared.Source;tol=max(1e-7,b.FrameDurationS*1e-4);index=0;
gray=zeros(prepared.Size);
if time<b.TimesS(1)-tol||time>=b.TimesS(end)+b.FrameDurationS-tol,return;end
index=find(b.TimesS<=time+tol,1,'last');
gray=Sample(prepared,double(b.Stack(:,:,index)),index);

function gray=Sample(p,imageData,index)
b=p.Source;x=p.X;y=p.Y;
if ~isempty(b.Drift)
    T=b.Drift(:,:,min(index,size(b.Drift,3)));inverse=T\eye(3);
    nx=x*inverse(1,1)+y*inverse(2,1)+inverse(3,1);
    y=x*inverse(1,2)+y*inverse(2,2)+inverse(3,2);x=nx;
end
gray=interp2(imageData,x,y,'linear',b.Limits(1));
gray=max(0,min(1,(gray-b.Limits(1))/(b.Limits(2)-b.Limits(1))));

function gray=Projection(p,autoContrast)
if nargin<2,autoContrast=false;end
native=max(p.Source.Stack,[],3);
if autoContrast
    limits=prctile(double(native(:)),[1 99.8]);
    if limits(2)<=limits(1),limits(2)=limits(1)+1;end
    p.Source.Limits=limits;
end
if isempty(p.Source.Drift)
    gray=Sample(p,double(native),1);
else
    gray=zeros(p.Size);
    for k=1:size(p.Source.Stack,3),gray=max(gray,Sample(p,double(p.Source.Stack(:,:,k)),k));end
end
