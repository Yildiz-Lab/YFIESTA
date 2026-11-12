function varargout=fTraceMovie(action,varargin)
%FTRACEMOVIE FIESTA molecule movies, no GUI/global-state dependencies.
% o=fTraceMovie('defaults');
% report=fTraceMovie('export',molecule,'trace.avi',o); % .avi/.mp4/.gif
% plan=fTraceMovie('prepare',molecule,o); state=fTraceMovie('state',plan,k);
% rgb=fTraceMovie('frame',plan,k); % the same renderer used by export
% Results columns: [frame, time in seconds, X in nm, Y in nm, ...].
switch lower(action)
    case 'defaults'
        varargout{1}=struct('FrameRate',30,'PlaybackSpeed',.25,'DecayMs',150, ...
            'BaselineAlpha',.20,'MarkerSize',2.5/3,'LineWidth',.5,'BufferFraction',.10, ...
            'CropNormalized',[],'Rotation','none','Theme','dark','Color',[], ...
            'TimeOrigin','relative','ScaleBarNm',[],'MaxDimension',960, ...
            'FrameSize',[],'ResolutionScale',1,'Quality',95,'Overwrite',false, ...
            'DeadTimeSpeedup',20,'ShowLegend',false,'ColorByChannel',false,'CropTiming','visible', ...
            'TimeRangeS',[],'Background',[],'ProgressFcn',[]);
    case 'prepare'
        varargout{1}=Prepare(varargin{:});
    case 'state'
        varargout{1}=FrameState(varargin{:});
    case 'frame'
        plan=varargin{1}; renderer=CreateRenderer(plan); cleanup=onCleanup(@()DisposeRenderer(renderer));
        [varargout{1:nargout}]=RenderFrame(renderer,plan,varargin{2});
    case 'renderer'
        varargout{1}=CreateRenderer(varargin{1});
    case 'render'
        [varargout{1:nargout}]=RenderFrame(varargin{:});
    case 'dispose'
        DisposeRenderer(varargin{1});
    case 'export'
        varargout{1}=Export(varargin{:});
    otherwise
        error('FIESTA:TraceMovieAction','Unknown trace movie action: %s',action);
end

function plan=Prepare(molecule,options)
if nargin<2,options=struct;end
o=fTraceMovie('defaults');
if ~isstruct(options)||~isscalar(options),error('FIESTA:TraceMovieOptions','Options must be a scalar struct.');end
fields=fieldnames(options);
for k=1:numel(fields)
    if ~isfield(o,fields{k}),error('FIESTA:TraceMovieOptions','Unknown option: %s',fields{k});end
    o.(fields{k})=options.(fields{k});
end
if ~isstruct(molecule)||isempty(molecule)||~isvector(molecule)||~isfield(molecule,'Results')
    error('FIESTA:TraceMovieMolecules','Pass a nonempty vector of FIESTA molecules.');
end
count=numel(molecule);traces=struct('Name',{},'Channel',{},'DataTimesS',{},'OriginalXYNm',{},'XYNm',{},'Color',{},'GlobalIndices',{});
motion=cell(count,1);offset=0;
for k=1:count
    r=double(molecule(k).Results);validateattributes(r,{'double'},{'2d','nonempty'});
    if size(r,2)<4,error('FIESTA:TraceMovieResults','Results needs time, X and Y columns (2, 3 and 4).');end
    if any(~isfinite(r(:,2:4)),'all')||any(diff(r(:,2))<0)
        error('FIESTA:TraceMovieResults','Time/X/Y must be finite, and timestamps must be nondecreasing (molecule %d).',k);
    end
    name=sprintf('Molecule %d',k);
    if isfield(molecule,'Name'),name=strtrim(regexprep(char(molecule(k).Name),'\s+',' '));end
    channel=1;if isfield(molecule,'Channel')&&~isempty(molecule(k).Channel),channel=molecule(k).Channel;end
    validateattributes(channel,{'numeric'},{'real','scalar','finite','integer','positive'});
    traces(k)=struct('Name',name,'Channel',double(channel),'DataTimesS',r(:,2),'OriginalXYNm',r(:,3:4), ...
        'XYNm',[],'Color',[],'GlobalIndices',offset+(1:size(r,1))');
    motion{k}=r(:,3:4)-mean(r(:,3:4),1);offset=offset+size(r,1);
end
for key={'FrameRate','PlaybackSpeed','DecayMs','MarkerSize','LineWidth','MaxDimension'}
    validateattributes(o.(key{1}),{'numeric'},{'real','scalar','finite','positive'});
end
validateattributes(o.FrameRate,{'numeric'},{'<=',120});
validateattributes(o.DeadTimeSpeedup,{'numeric'},{'real','scalar','finite','>=',1});
validateattributes(o.BaselineAlpha,{'numeric'},{'real','scalar','finite','>=',0,'<=',1});
validateattributes(o.BufferFraction,{'numeric'},{'real','scalar','finite','nonnegative'});
validateattributes(o.Quality,{'numeric'},{'real','scalar','finite','>=',0,'<=',100});
validateattributes(o.MaxDimension,{'numeric'},{'integer','>=',320});
validateattributes(o.ResolutionScale,{'numeric'},{'real','scalar','finite','integer','>=',1,'<=',4});
o.Rotation=validatestring(o.Rotation,{'none','horizontal','vertical'});
o.Theme=validatestring(o.Theme,{'dark','light'});
o.TimeOrigin=validatestring(o.TimeOrigin,{'relative','absolute'});
o.CropTiming=validatestring(o.CropTiming,{'visible','full'});
if ~isempty(o.ProgressFcn)&&~isa(o.ProgressFcn,'function_handle')
    error('FIESTA:TraceMovieProgress','ProgressFcn must be a function handle.');
end
if ~isempty(o.Background),o.Background=fTraceMovieBackground('validate',o.Background);end
validateattributes(o.Overwrite,{'logical','numeric'},{'scalar','binary'});
validateattributes(o.ShowLegend,{'logical','numeric'},{'scalar','binary'});
validateattributes(o.ColorByChannel,{'logical','numeric'},{'scalar','binary'});
colors=AssignColors(traces,o);
legendLabels={traces.Name};legendColors=colors;
if o.ColorByChannel&&isempty(o.Color)
    [channels,first]=unique([traces.Channel]);legendColors=colors(first,:);
    legendLabels=arrayfun(@(c)sprintf('Channel %g',c),channels,'UniformOutput',false);
end
xy=vertcat(traces.OriginalXYNm);center=mean(xy,1);R=eye(2);motion=vertcat(motion{:});
if ~strcmp(o.Rotation,'none') && any(motion~=0,'all')
    % Within-trace motion determines the axis, rather than distance between traces.
    [~,~,v]=svd(motion,0);axisVector=v(:,1);direction=[0 0];
    for k=1:count,direction=direction+traces(k).OriginalXYNm(end,:)-traces(k).OriginalXYNm(1,:);end
    if dot(axisVector,direction)<0 || (norm(direction)==0 && axisVector(1)<0),axisVector=-axisVector;end
    perpendicular=[-axisVector(2);axisVector(1)];
    if strcmp(o.Rotation,'horizontal'),R=[axisVector perpendicular];else,R=[-perpendicular axisVector];end
end
xy=(xy-center)*R;
for k=1:count,traces(k).XYNm=xy(traces(k).GlobalIndices,:);traces(k).Color=colors(k,:);end
low=min(xy,[],1);high=max(xy,[],1);minimumSpan=50;
% Keep enough image context to crop a rotated, almost straight trace.
if ~isempty(o.Background),minimumSpan=max(minimumSpan,20*o.Background.PixelSizeNm);end
span=max(high-low,minimumSpan);mid=(low+high)/2;
fullBounds=[mid-span*(.5+o.BufferFraction);mid+span*(.5+o.BufferFraction)];
bounds=fullBounds;
if ~isempty(o.CropNormalized)
    box=o.CropNormalized;
    validateattributes(box,{'numeric'},{'real','finite','vector','numel',4,'>=',0,'<=',1});
    box=reshape(box,1,4);
    if any(box(3:4)<=0)||box(1)+box(3)>1+eps||box(2)+box(4)>1+eps
        error('FIESTA:TraceMovieCrop','Crop is [left top width height], within [0,1], with positive size.');
    end
    fullSpan=fullBounds(2,:)-fullBounds(1,:);
    bounds=[fullBounds(1,:)+box(1:2).*fullSpan;fullBounds(1,:)+(box(1:2)+box(3:4)).*fullSpan];
end
% Image coordinates: X rightwards, Y downwards, as in FIESTA's image view.
span=bounds(2,:)-bounds(1,:);margin=24;header=42;footer=50;
if isempty(o.FrameSize)
    scale=(o.MaxDimension-2*margin-header-footer)/max(span);
    plotSize=max(160,ceil(span*scale-1e-9));frameSize=2*ceil((plotSize+[2*margin 2*margin+header+footer])/2);
    frameSize(1)=max(320,frameSize(1)); % room for the title and millisecond clock
else
    validateattributes(o.FrameSize,{'numeric'},{'real','finite','integer','vector','numel',2,'>=',240});
    frameSize=2*ceil(reshape(o.FrameSize,1,2)/2);
    if frameSize(1)<320,error('FIESTA:TraceMovieFrameSize','Frame width must be at least 320 pixels for the title and clock.');end
end
legendColumns=max(1,min(4,floor((frameSize(1)-2*margin)/120)));legendHeight=0;
if o.ShowLegend&&(count>1||o.ColorByChannel),legendHeight=20*ceil(numel(legendLabels)/legendColumns);end
if isempty(o.FrameSize),frameSize(2)=frameSize(2)+legendHeight;end
% macOS MPEG-4 encoding can corrupt RGB row strides at non-block-aligned sizes.
% Use the same 16-pixel-aligned canvas for preview and all export formats.
frameSize=16*ceil(frameSize/16);
footer=footer+legendHeight;
available=frameSize-[2*margin 2*margin+header+footer];
if any(available<40),error('FIESTA:TraceMovieFrameSize','Frame is too small for the legend. Increase FrameSize or set ShowLegend=false.');end
scale=min(available./span);plotSize=span*scale;
rect=[margin+(available(1)-plotSize(1))/2 footer+margin+(available(2)-plotSize(2))/2 plotSize];
bar=o.ScaleBarNm;
if isempty(bar)
    footerWidth=frameSize(1)-2*margin;
    target=min(.35*footerWidth,max(80,.22*footerWidth))/scale;
    power=10^floor(log10(target));choices=[1 2 5]*power;
    bar=max(choices(choices<=target));
end
validateattributes(bar,{'numeric'},{'real','scalar','finite','positive'});
if bar*scale>.8*(frameSize(1)-2*margin),error('FIESTA:TraceMovieScaleBar','Scale bar must fit within the movie footer.');end
t=vertcat(traces.DataTimesS);timeOrigin=min(t);timeRange=[min(t) max(t)];
if ~isempty(o.TimeRangeS)
    validateattributes(o.TimeRangeS,{'numeric'},{'real','finite','vector','numel',2});
    if o.TimeRangeS(1)>o.TimeRangeS(2),error('FIESTA:TraceMovieTime','Time range must be increasing.');end
    timeRange=[max(timeRange(1),o.TimeRangeS(1)) min(timeRange(2),o.TimeRangeS(2))];
end
if ~isempty(o.CropNormalized)&&strcmp(o.CropTiming,'visible')
    inside=all(xy>=bounds(1,:) & xy<=bounds(2,:),2)&t>=timeRange(1)&t<=timeRange(2);
    if ~any(inside),error('FIESTA:TraceMovieEmptyCrop','The crop contains no recorded localizations in the chosen time range.');end
    timeRange=[min(t(inside)) max(t(inside))];
end
if timeRange(1)>timeRange(2),error('FIESTA:TraceMovieTime','Time range does not overlap the selected molecules.');end
[frameTimes,timeline,movieTimes]=BuildTimeline(traces,o,timeRange);
legendStart=-inf(1,numel(legendLabels));
if ~o.ColorByChannel||~isempty(o.Color),legendStart=TraceAppearanceTimes(traces,bounds);end
% Scale the complete layout; keep physical bounds, clock and legend slots fixed.
frameSize=frameSize*o.ResolutionScale;rect=rect*o.ResolutionScale;
scale=scale*o.ResolutionScale;legendHeight=legendHeight*o.ResolutionScale;
duration=diff(timeRange);n=numel(frameTimes);name=traces(1).Name;
if count>1,name=sprintf('%d molecules',count);end
plan=struct('Name',name,'Options',o,'DataTimesS',t,'OriginalXYNm',vertcat(traces.OriginalXYNm), ...
    'XYNm',xy,'CenterNm',center,'RotationMatrix',R,'FullBoundsNm',fullBounds, ...
    'BoundsNm',bounds,'FrameTimesS',frameTimes,'FrameCount',n,'FrameSize',frameSize, ...
    'PlotRectPx',rect,'PixelsPerNm',scale,'ScaleBarNm',bar,'SourceDurationS',duration, ...
    'Traces',traces,'MoleculeCount',count,'Colors',colors,'Timeline',timeline, ...
    'LegendColumns',legendColumns,'LegendHeight',legendHeight,'LegendLabels',{legendLabels}, ...
    'LegendColors',legendColors,'LegendStartTimesS',legendStart,'TimeOriginS',timeOrigin, ...
    'TimeRangeS',timeRange,'MovieTimesS',movieTimes);

function times=TraceAppearanceTimes(traces,bounds)
times=inf(1,numel(traces));
for k=1:numel(traces)
    xy=traces(k).XYNm;t=traces(k).DataTimesS;
    first=find(all(xy>=bounds(1,:) & xy<=bounds(2,:),2),1);
    if ~isempty(first),times(k)=t(first);end
    % A connecting segment can cross the crop with both endpoints outside it.
    % The renderer reveals that segment only after its second localization.
    [lo,hi]=ClipSegments(xy(1:end-1,:),diff(xy,1,1),bounds(1,:),bounds(2,:));
    first=find(lo<=hi,1);if ~isempty(first),times(k)=min(times(k),t(first+1));end
end

function [lo,hi]=ClipSegments(a,v,low,high)
lo=zeros(size(a,1),1);hi=ones(size(a,1),1);
for dim=1:2
    moving=v(:,dim)~=0;
    ends=([low(dim) high(dim)]-a(moving,dim))./v(moving,dim);
    lo(moving)=max(lo(moving),min(ends,[],2));hi(moving)=min(hi(moving),max(ends,[],2));
    hi(~moving & (a(:,dim)<low(dim)|a(:,dim)>high(dim)))=-1;
end

function colors=AssignColors(traces,o)
count=numel(traces);
if ~isempty(o.Color)
    colors=o.Color;
    validateattributes(colors,{'numeric'},{'real','finite','2d','>=',0,'<=',1});
    if isvector(colors)&&numel(colors)==3,colors=reshape(colors,1,3);end
    if size(colors,2)~=3||~ismember(size(colors,1),[1 count])
        error('FIESTA:TraceMovieColor','Color must be one RGB triplet or one RGB row per molecule.');
    end
    if size(colors,1)==1,colors=repmat(colors,count,1);end
    return;
end
% Okabe-Ito hues, reused with local spatial/time contrast rather than global ties.
palette=[86 180 233;230 159 0;0 158 115;204 121 167;213 94 0;0 114 178;240 228 66]/255;
if strcmp(o.Theme,'light')
    palette=palette([6 5 3 4 2 1 7],:);
    palette(end,:)=palette(end,:)*.65; % yellow needs more contrast on white
end
if o.ColorByChannel
    colors=palette(1+mod([traces.Channel]-1,size(palette,1)),:);return;
end
if count==1,colors=palette(1,:);return;end
[~,canonical]=sort(string({traces.Name}));traces=traces(canonical);
samples=cell(count,1);low=zeros(count,2);high=low;center=low;interval=zeros(count,2);
for k=1:count
    xy=traces(k).OriginalXYNm;
    samples{k}=xy(unique(round(linspace(1,size(xy,1),min(64,size(xy,1))))),:);
    low(k,:)=min(xy,[],1);high(k,:)=max(xy,[],1);center(k,:)=mean(xy,1);
    interval(k,:)=traces(k).DataTimesS([1 end]);
end
scale=max(1,median(vecnorm(high-low,2,2)));
timeScale=max(3*o.DecayMs/1000,median(diff(interval,1,2)));
neighbors=min(24,count-1);row=zeros(count*neighbors,1);column=row;values=row;
for k=1:count
    gap=max(0,max(interval(:,1),interval(k,1))-min(interval(:,2),interval(k,2)));
    temporal=o.BaselineAlpha+(1-o.BaselineAlpha)*exp(-gap/timeScale);
    boxDistance=vecnorm(max(0,max(low-high(k,:),low(k,:)-high)),2,2);
    score=temporal./(1+(boxDistance/scale).^2+.01*sum((center-center(k,:)).^2,2)/scale^2);
    score(k)=-Inf;[~,candidate]=maxk(score,neighbors);
    for n=1:neighbors
        j=candidate(n);
        a=samples{k};b=samples{j};
        d=hypot(a(:,1)-b(:,1)',a(:,2)-b(:,2)');
        idx=(k-1)*neighbors+n;row(idx)=k;column(idx)=j;
        values(idx)=temporal(j)/(1+(min(d,[],'all')/scale)^2);
    end
end
weight=sparse(row,column,values,count,count);weight=max(weight,weight');
lab=rgb2lab(reshape(palette,[],1,3));lab=reshape(lab,[],3);
penalty=exp(-sum((reshape(lab,[],1,3)-reshape(lab,1,[],3)).^2,3)/(2*35^2));
assigned=zeros(count,1);usage=zeros(size(palette,1),1);priority=full(sum(weight,2));
for step=1:count
    [~,k]=max(priority);[~,choice]=min(ColorCosts(k,weight,assigned,usage,penalty));
    assigned(k)=choice;usage(choice)=usage(choice)+1;priority(k)=-Inf;
    [~,near,w]=find(weight(k,:));open=assigned(near)==0;
    priority(near(open))=priority(near(open))+reshape(w(open),[],1);
end
% Refine against all local neighbors after the first greedy pass.
for pass=1:2
    for k=1:count
        old=assigned(k);usage(old)=usage(old)-1;
        scores=ColorCosts(k,weight,assigned,usage,penalty);[best,choice]=min(scores);
        if scores(old)<=best+1e-12,choice=old;end
        assigned(k)=choice;usage(choice)=usage(choice)+1;
    end
end
colors=zeros(count,3);colors(canonical,:)=palette(assigned,:);

function scores=ColorCosts(k,weight,assigned,usage,penalty)
[~,near,w]=find(weight(k,:));keep=assigned(near)>0;near=near(keep);w=w(keep);
scores=zeros(numel(usage),1);
if ~isempty(near),scores=penalty(:,assigned(near))*w(:);end
scores=scores+.02*max(1,sum(w))*usage/max(1,numel(assigned)/numel(usage));
if numel(assigned)<=numel(usage),scores(usage>0)=Inf;end

function [times,timeline,movieTimes]=BuildTimeline(traces,o,timeRange)
intervals=zeros(numel(traces),2);
for k=1:numel(traces),intervals(k,:)=traces(k).DataTimesS([1 end]);end
intervals=sortrows(intervals,1);merged=intervals(1,:);
for k=2:size(intervals,1)
    if intervals(k,1)<=merged(end,2),merged(end,2)=max(merged(end,2),intervals(k,2));
    else,merged(end+1,:)=intervals(k,:);end %#ok<AGROW>
end
edges=unique([merged(:);timeRange(:)]);edges=edges(edges>=timeRange(1)&edges<=timeRange(2));
timeline=struct('SourceStartS',{},'SourceEndS',{},'IsDeadTime',{},'SpeedFactor',{}, ...
    'MovieStartS',{},'MovieEndS',{},'FirstFrame',{},'LastFrame',{});
movieEdges=zeros(size(edges));
for k=1:numel(edges)-1
    a=edges(k);b=edges(k+1);mid=(a+b)/2;
    idle=~any(mid>=intervals(:,1)&mid<=intervals(:,2));factor=1;
    if idle,factor=o.DeadTimeSpeedup;end
    movieEdges(k+1)=movieEdges(k)+(b-a)/(o.PlaybackSpeed*factor);
    timeline(k)=struct('SourceStartS',a,'SourceEndS',b,'IsDeadTime',idle, ...
        'SpeedFactor',factor,'MovieStartS',movieEdges(k),'MovieEndS',movieEdges(k+1), ...
        'FirstFrame',0,'LastFrame',0);
end
% One continuous physical-time -> playback-time map, sampled at output FPS.
% Molecule boundaries and irregular localization spacing add no output time.
n=max(1,ceil(movieEdges(end)*o.FrameRate)+1);
if ~isfinite(n)||n>1e6,error('FIESTA:TraceMovieLength','Export would exceed one million output images. Increase PlaybackSpeed.');end
movieTimes=min((0:n-1)'/o.FrameRate,movieEdges(end));
if isscalar(edges),times=edges;
else
    [uniqueMovie,indices]=unique(movieEdges,'last');
    times=interp1(uniqueMovie,edges(indices),movieTimes,'linear');
    times(1)=timeRange(1);times(end)=timeRange(2);
end
for k=1:numel(timeline)
    indices=find(movieTimes>=timeline(k).MovieStartS & movieTimes<=timeline(k).MovieEndS);
    if ~isempty(indices),timeline(k).FirstFrame=indices(1);timeline(k).LastFrame=indices(end);end
end

function state=FrameState(plan,k)
validateattributes(k,{'numeric'},{'real','finite','scalar','integer','>=',1,'<=',plan.FrameCount});
t=plan.FrameTimesS(k);indices=find(plan.DataTimesS<=t);
ageMs=max(0,(t-plan.DataTimesS(indices))*1000);
a=plan.Options.BaselineAlpha+(1-plan.Options.BaselineAlpha)*exp(-ageMs/plan.Options.DecayMs);
label=t;if strcmp(plan.Options.TimeOrigin,'relative'),label=t-plan.TimeOriginS;end
traceStates=struct('Name',{},'VisibleIndices',{},'AgeMs',{},'Alpha',{},'IsActive',{});active=false;
for j=1:plan.MoleculeCount
    trace=plan.Traces(j);local=find(trace.DataTimesS<=t);ages=max(0,(t-trace.DataTimesS(local))*1000);
    isActive=t>=trace.DataTimesS(1)&&t<=trace.DataTimesS(end);active=active||isActive;
    traceStates(j)=struct('Name',trace.Name,'VisibleIndices',local,'AgeMs',ages, ...
        'Alpha',plan.Options.BaselineAlpha+(1-plan.Options.BaselineAlpha)*exp(-ages/plan.Options.DecayMs), ...
        'IsActive',isActive);
end
state=struct('TimeS',t,'LabelTimeS',label,'VisibleIndices',indices,'AgeMs',ageMs,'Alpha',a, ...
    'Traces',traceStates,'IsDeadTime',~active,'SpeedFactor',1,'LegendVisible',t>=plan.LegendStartTimesS);
if ~active,state.SpeedFactor=plan.Options.DeadTimeSpeedup;end

function renderer=CreateRenderer(plan)
o=plan.Options;if strcmp(o.Theme,'dark'),bg=[.045 .055 .075];fg=[.94 .95 .97];else,bg=[1 1 1];fg=[.08 .09 .11];end
old=get(0,'CurrentFigure');
% Keep the hidden figure at base size, so high-resolution exports are not
% constrained by the monitor. Raster layers and print DPI supply the extra pixels.
s=o.ResolutionScale;baseSize=plan.FrameSize/s;dpi=round(get(0,'ScreenPixelsPerInch'));
f=figure('Visible','off','Units','pixels','Position',[40 40 baseSize], ...
    'Color',bg,'MenuBar','none','ToolBar','none','NumberTitle','off','Name','Trace movie render', ...
    'Resize','off','HandleVisibility','off','InvertHardcopy','off');
resolution='-r0';
if s>1
    set(f,'PaperUnits','inches','PaperPositionMode','manual', ...
        'PaperPosition',[0 0 baseSize/dpi],'PaperSize',baseSize/dpi);
    resolution=sprintf('-r%d',dpi*s);
end
try
    ax=axes('Parent',f,'Units','pixels','Position',plan.PlotRectPx/s,'Color',bg, ...
        'XLim',plan.BoundsNm(:,1)','YLim',plan.BoundsNm(:,2)','YDir','reverse', ...
        'DataAspectRatio',[1 1 1],'Visible','off','NextPlot','add');
    raster=cell(1,plan.MoleculeCount);
    for j=1:plan.MoleculeCount,raster{j}=PrepareRaster(plan,plan.Traces(j));end
    rasterSize=raster{1}.Size;color=repmat(reshape(plan.Colors(1,:),1,1,3),rasterSize);
    span=diff(plan.BoundsNm,1,1);pixel=span./rasterSize([2 1]);
    background=[];
    if ~isempty(o.Background)
        background=fTraceMovieBackground('prepare',o.Background,plan,rasterSize);
    end
    layer=image(ax,'XData',[plan.BoundsNm(1,1)+pixel(1)/2 plan.BoundsNm(2,1)-pixel(1)/2], ...
        'YData',[plan.BoundsNm(1,2)+pixel(2)/2 plan.BoundsNm(2,2)-pixel(2)/2], ...
        'CData',color);
    W=baseSize(1);H=baseSize(2);
    annotation(f,'textbox',[24/W (H-36)/H .56 28/H],'String',plan.Name,'Interpreter','none', ...
        'EdgeColor','none','Color',fg,'FontSize',13,'FontName','Helvetica','VerticalAlignment','middle');
    clock=annotation(f,'textbox',[.59 (H-36)/H .38 28/H],'String','','Interpreter','none', ...
        'EdgeColor','none','Color',fg,'FontSize',14,'FontName','Helvetica','HorizontalAlignment','right');
    lengthPx=plan.ScaleBarNm*plan.PixelsPerNm/s;right=W-24;
    annotation(f,'line',[right-lengthPx right]/W,[34 34]/H,'Color',fg,'LineWidth',2);
    if plan.ScaleBarNm>=1000,label=sprintf('%g %sm',plan.ScaleBarNm/1000,char(181));else,label=sprintf('%g nm',plan.ScaleBarNm);end
    labelWidth=max(100,lengthPx);
    annotation(f,'textbox',[(right-lengthPx/2-labelWidth/2)/W 6/H labelWidth/W 25/H],'String',label,'Interpreter','none', ...
        'EdgeColor','none','Color',fg,'FontSize',11,'HorizontalAlignment','center','Margin',0);
    speedLabel=annotation(f,'textbox',[24/W 25/H .5 25/H],'String','','Interpreter','none', ...
        'EdgeColor','none','Color',fg,'FontSize',10,'Margin',0);
    legend=gobjects(2,0);
    if plan.LegendHeight>0
        columnWidth=(W-48)/plan.LegendColumns;
        legend=gobjects(2,numel(plan.LegendLabels));
        for j=1:numel(plan.LegendLabels)
            x=24+mod(j-1,plan.LegendColumns)*columnWidth;y=52+plan.LegendHeight/s-20*ceil(j/plan.LegendColumns);
            legend(1,j)=annotation(f,'line',[x x+15]/W,[y+8 y+8]/H,'Color',plan.LegendColors(j,:),'LineWidth',2,'Visible','off');
            label=regexprep(plan.LegendLabels{j},'^Molecule\s+','');
            legend(2,j)=annotation(f,'textbox',[(x+20)/W y/H (columnWidth-20)/W 20/H],'String',label, ...
                'Interpreter','none','EdgeColor','none','Color',fg,'FontSize',10,'Margin',0,'Visible','off');
        end
    end
    renderer=struct('Figure',f,'PreviousFigure',old,'Layer',layer,'Raster',{raster}, ...
        'Clock',clock,'SpeedLabel',speedLabel,'Background',bg,'Foreground',fg, ...
        'StackBackground',background,'Legend',legend,'PrintResolution',resolution);
catch err
    delete(f);rethrow(err);
end

function [rgb,state]=RenderFrame(renderer,plan,k)
state=FrameState(plan,k);sizePx=renderer.Raster{1}.Size;
state.BackgroundFrame=0;
baseImage=repmat(reshape(renderer.Background,1,1,3),sizePx);
if ~isempty(renderer.StackBackground)
    [imageData,state.BackgroundFrame]=fTraceMovieBackground('frame',renderer.StackBackground,state.TimeS);
    baseImage=repmat(imageData,1,1,3);
end
alpha=zeros(sizePx);colorSum=zeros([sizePx 3]);ties=zeros(sizePx);
if nargout>1,state.TraceCoverageAlpha=cell(1,plan.MoleculeCount);end
for j=1:plan.MoleculeCount
    r=renderer.Raster{j};
    points=Paint(r.Points,state.TimeS,plan.Options,r.Size,r.Supersample,plan.Options.MarkerSize*plan.Options.ResolutionScale/2);
    lines=Paint(r.Lines,state.TimeS,plan.Options,r.Size,r.Supersample,plan.Options.LineWidth*plan.Options.ResolutionScale/2);
    % Max-composite within each trace prevents dense old points building opacity.
    a=imresize(max(points,lines),r.Size,'box');
    if nargout>1,state.TraceCoverageAlpha{j}=a;end
    wins=a>alpha;equal=a==alpha & a>0;
    for c=1:3
        channel=colorSum(:,:,c);channel(wins)=plan.Colors(j,c);
        channel(equal)=channel(equal)+plan.Colors(j,c);colorSum(:,:,c)=channel;
    end
    ties(wins)=1;ties(equal)=ties(equal)+1;alpha=max(alpha,a);
end
% At crossings the brighter/recent trace wins; exact ties blend symmetrically.
state.CoverageAlpha=alpha;
% One opaque image avoids depth/transparency artifacts across native ROI windows.
set(renderer.Layer,'CData',(colorSum./max(ties,1)).*alpha+baseImage.*(1-alpha));
set(renderer.Clock,'String',sprintf('t = %.3f s',state.LabelTimeS));
speed='';if state.IsDeadTime&&state.SpeedFactor>1,speed=sprintf('Idle %gx',state.SpeedFactor);end
set(renderer.SpeedLabel,'String',speed);
if ~isempty(renderer.Legend)
    set(renderer.Legend(:,state.LegendVisible),'Visible','on');
    set(renderer.Legend(:,~state.LegendVisible),'Visible','off');
end
drawnow;
% Scene-based capture keeps pixels consistent with a native crop editor open.
rgb=print(renderer.Figure,'-RGBImage',renderer.PrintResolution,'-painters');
if ~isequal([size(rgb,2) size(rgb,1)],plan.FrameSize),rgb=imresize(rgb,plan.FrameSize([2 1]),'bicubic');end

function raster=PrepareRaster(plan,trace)
ss=2;sizePx=max(1,round(plan.PlotRectPx([4 3])));grid=sizePx*ss;
span=diff(plan.BoundsNm,1,1);
xy=(trace.XYNm-plan.BoundsNm(1,:))./span.*grid([2 1])+.5;
points=[xy trace.DataTimesS trace.DataTimesS];
segments=cell(max(0,size(xy,1)-1),1);
[lo,hi]=ClipSegments(xy(1:end-1,:),diff(xy,1,1),[.5 .5],grid([2 1])+.5);
for k=1:numel(segments)
    a=xy(k,:);b=xy(k+1,:);v=b-a;
    % Clip before sampling so very tight crops do not create huge arrays.
    if lo(k)>hi(k),continue;end
    fractions=linspace(lo(k),hi(k),max(2,ceil(norm(v)*(hi(k)-lo(k))*2)+1))';
    times=trace.DataTimesS(k)+fractions*(trace.DataTimesS(k+1)-trace.DataTimesS(k));
    segments{k}=[a+fractions.*v times repmat(trace.DataTimesS(k+1),numel(times),1)];
end
if isempty(segments),lines=zeros(0,4);else,lines=vertcat(segments{:});end
if isempty(lines),lines=zeros(0,4);end % every segment may have been clipped out
raster=struct('Size',sizePx,'Supersample',ss,'Points',points,'Lines',lines);

function mask=Paint(samples,time,options,sizePx,ss,radius)
grid=sizePx*ss;radius=radius*ss;
samples=samples(samples(:,4)<=time,:);
pad=ceil(radius+1);
samples=samples(samples(:,1)>=-pad & samples(:,1)<=grid(2)+pad & ...
    samples(:,2)>=-pad & samples(:,2)<=grid(1)+pad,:);
if isempty(samples),mask=zeros(grid);return;end
alpha=options.BaselineAlpha+(1-options.BaselineAlpha)*exp(-max(0,(time-samples(:,3))*1000)/options.DecayMs);
[dx,dy]=meshgrid(-pad:pad,-pad:pad);
x=round(samples(:,1))'+dx(:);y=round(samples(:,2))'+dy(:);
coverage=max(0,min(1,radius+.5-hypot(x-samples(:,1)',y-samples(:,2)'))).*alpha';
valid=x>=1 & x<=grid(2) & y>=1 & y<=grid(1) & coverage>0;
indices=y(valid)+(x(valid)-1)*grid(1);
mask=reshape(accumarray(indices,coverage(valid),[prod(grid) 1],@max,0),grid);

function DisposeRenderer(renderer)
if ishghandle(renderer.Figure),delete(renderer.Figure);end
if ~isempty(renderer.PreviousFigure)&&ishghandle(renderer.PreviousFigure),set(0,'CurrentFigure',renderer.PreviousFigure);end

function report=Export(molecule,file,options)
if nargin<3,options=struct;end
plan=Prepare(molecule,options);o=plan.Options;file=char(file);
[folder,~,ext]=fileparts(file);if isempty(folder),folder=pwd;file=fullfile(folder,file);end
if ~isfolder(folder),error('FIESTA:TraceMovieFolder','Output folder does not exist: %s',folder);end
if isfile(file)&&~o.Overwrite,error('FIESTA:TraceMovieExists','Output exists; choose another name or set Overwrite=true.');end
ext=lower(ext);if ~ismember(ext,{'.avi','.mp4','.gif'}),error('FIESTA:TraceMovieFormat','Choose .avi, .mp4 or .gif.');end
if strcmp(ext,'.gif')&&o.FrameRate>50,error('FIESTA:TraceMovieGifRate','GIF export supports up to 50 fps (centisecond timing).');end
if strcmp(ext,'.mp4')
    profiles=VideoWriter.getProfiles;
    if ~any(strcmp({profiles.Name},'MPEG-4')),error('FIESTA:TraceMovieEncoder','MPEG-4 is unavailable in this MATLAB installation. Use AVI or GIF.');end
end
temporary=[tempname(folder) ext];tempCleanup=onCleanup(@()DeleteTemporary(temporary));
renderer=CreateRenderer(plan);renderCleanup=onCleanup(@()DisposeRenderer(renderer));
writer=[];
try
    if strcmp(ext,'.gif')
        base=unique(plan.Colors,'rows','stable');
        if size(base,1)>32,base=base(round(linspace(1,size(base,1),32)),:);end
        shadeCount=128;textCount=128;
        if size(base,1)>1,shadeCount=floor(224/size(base,1));textCount=32;end
        if ~isempty(o.Background),shadeCount=floor(128/size(base,1));textCount=128;end
        levels=linspace(0,1,shadeCount)';palette=[];
        for j=1:size(base,1)
            palette=[palette;renderer.Background+(base(j,:)-renderer.Background).*levels]; %#ok<AGROW>
        end
        gray=renderer.Background+(renderer.Foreground-renderer.Background).*linspace(0,1,textCount)';
        if ~isempty(o.Background),gray=repmat(linspace(0,1,textCount)',1,3);end
        palette=[palette;gray];
        palette=unique(palette,'rows','stable');
        edges=round((0:plan.FrameCount)*100/o.FrameRate);delays=diff(edges)/100;
    else
        profile='Motion JPEG AVI';if strcmp(ext,'.mp4'),profile='MPEG-4';end
        writer=VideoWriter(temporary,profile);writer.FrameRate=o.FrameRate;writer.Quality=o.Quality;open(writer);
    end
    for k=1:plan.FrameCount
        if ~isempty(o.ProgressFcn)&&~o.ProgressFcn(k-1,plan.FrameCount)
            error('FIESTA:TraceMovieCancelled','Export cancelled.');
        end
        rgb=RenderFrame(renderer,plan,k);
        if strcmp(ext,'.gif')
            indexed=rgb2ind(rgb,palette,'nodither');
            if k==1
                imwrite(indexed,palette,temporary,'gif','LoopCount',Inf,'DelayTime',delays(k));
            else
                imwrite(indexed,palette,temporary,'gif','WriteMode','append','DelayTime',delays(k));
            end
        else
            writeVideo(writer,rgb);
        end
        if mod(k,max(1,ceil(plan.FrameCount/10)))==0 || k==plan.FrameCount
            fprintf('Trace movie %s: %d/%d frames\n',ext,k,plan.FrameCount);
        end
    end
    if ~isempty(writer),close(writer);writer=[];end
    if ~isempty(o.ProgressFcn)&&~o.ProgressFcn(plan.FrameCount,plan.FrameCount)
        error('FIESTA:TraceMovieCancelled','Export cancelled.');
    end
    if isfile(file)&&~o.Overwrite,error('FIESTA:TraceMovieExists','Output appeared during export; it has not been replaced.');end
    [ok,message]=movefile(temporary,file,'f');if ~ok,error('FIESTA:TraceMovieSave','%s',message);end
catch err
    if ~isempty(writer)
        try
            close(writer);
        catch
        end
    end
    rethrow(err);
end
report=struct('File',file,'MoleculeName',plan.Name,'FrameCount',plan.FrameCount, ...
    'FrameSize',plan.FrameSize,'FrameRate',o.FrameRate,'PlaybackSpeed',o.PlaybackSpeed, ...
    'SourceDurationS',plan.SourceDurationS,'MovieDurationS',plan.FrameCount/o.FrameRate, ...
    'ScaleBarNm',plan.ScaleBarNm,'Options',o,'MoleculeCount',plan.MoleculeCount, ...
    'MoleculeNames',{{plan.Traces.Name}},'Channels',[plan.Traces.Channel],'Colors',plan.Colors,'Timeline',plan.Timeline, ...
    'DeadTimeSpeedup',o.DeadTimeSpeedup);

function DeleteTemporary(file)
if isfile(file),delete(file);end
