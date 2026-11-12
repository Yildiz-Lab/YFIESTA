function varargout=fTraceMovieGui(action,varargin)
%FTRACEMOVIEGUI One-window selected-molecule exporter; data are copied at open.
% create(molecules,backgroundCell) also supports isolated validation without globals.
switch lower(action)
    case ''
    case 'create'
        h=Create(varargin{:});if nargout,varargout{1}=h;end
    case 'install'
        main=getappdata(0,'hMainGui');
        if isempty(main)||~isgraphics(main.fig),return;end
        item=findall(main.Menu.mView,'Tag','mExportSelectedMovie');
        if isempty(item)
            item=uimenu(main.Menu.mView,'Label','Export selected molecule movie', ...
                'Tag','mExportSelectedMovie','Separator','on','Callback',@(s,e)fTraceMovieGui('create'));
        end
        main.Menu.mExportSelectedMovie=item;setappdata(0,'hMainGui',main);
    case 'options',varargout{1}=ReadOptions(varargin{1});
    case 'context',[varargout{1:nargout}]=CaptureContext;
    case 'refresh',Refresh(varargin{1});
    case 'export',Export(varargin{1});
    case 'crop',StoreCrop(varargin{1});
end

function h=Create(varargin)
if isempty(varargin)
    [m,sources,folder]=CaptureContext;
else
    m=varargin{1};sources={};folder=pwd;
    if numel(varargin)>1,sources=varargin{2};end
    if isstruct(sources),sources={sources};end
end
h=[];
if isempty(m)
    errordlg('Select one or more molecules in FIESTA, then choose this command again.','Export selected molecule movie');return;
end
old=findall(0,'Tag','TraceMovieDialog');
for k=1:numel(old)
    state=getappdata(old(k),'UI');
    if ~isempty(state)&&state.Busy,figure(old(k));return;end
    Close(old(k));
end
remembered=getappdata(0,'TraceMovieOutputFolder');if ischar(remembered)&&isfolder(remembered),folder=remembered;end
c=[50 50 50]/255;fg=[.9 .9 .9];editColor=[.12 .12 .12];button=[.16 .16 .16];
h=figure('Name','Export selected molecule movie','Tag','TraceMovieDialog','Units','pixels', ...
    'Position',[100 80 1180 800],'MenuBar','none','ToolBar','none','NumberTitle','off', ...
    'Color',c,'CloseRequestFcn',@(s,e)Close(s),'DefaultUicontrolBackgroundColor',c, ...
    'DefaultUicontrolForegroundColor',fg,'Visible','off','Resize','on');movegui(h,'center');
ui=struct('Molecules',m,'Sources',{sources},'Busy',false,'Cancel',false,'Closing',false, ...
    'Playing',false,'Dirty',true,'Plan',[],'BasePlan',[],'ROI',[],'Renderer',[],'Rotation','none');
uicontrol(h,'Style','text','Units','normalized','Position',[.02 .945 .96 .04], ...
    'String',sprintf('%d selected molecule(s) — selection captured when this window opened',numel(m)), ...
    'HorizontalAlignment','left','FontSize',12);
settings=uitabgroup(h,'Units','normalized','Position',[.02 .155 .35 .78]);
traceTab=uitab(settings,'Title','Traces','BackgroundColor',c);
backgroundTab=uitab(settings,'Title','Background & crop','BackgroundColor',c);
exportTab=uitab(settings,'Title','Export','BackgroundColor',c);
o=fTraceMovie('defaults');
keys={'PlaybackSpeed','FrameRate','DecayMs','BaselineAlpha','MarkerSize','LineWidth','DeadTimeSpeedup','BufferFraction','MaxDimension','ScaleBarNm'};
labels={'Source seconds / movie second','Movie FPS (output sampling)','Fade decay (ms)','Old-trace opacity (0–1)','Point diameter (base pixels)', ...
    'Connecting line width (base pixels)','Dead-time speedup (1 = normal)','Padding fraction','Base canvas size (pixels)','Scale bar (nm; blank = auto)'};
for k=1:numel(keys)
    key=keys{k};y=.92-(k-1)*.071;
    Label(traceTab,labels{k},[.04 y .68 .047]);
    ui.(key)=Edit(traceTab,num2str(o.(key)),[.74 y .22 .052],key,@(s,e)Dirty(h));
end
Label(traceTab,'Dominant motion axis',[.04 .205 .55 .045]);
ui.Rotation=Popup(traceTab,{'No rotation','Horizontal','Vertical'},[.61 .205 .35 .052],'Rotation',@(s,e)RotationChanged(h));
Label(traceTab,'Movie theme',[.04 .145 .55 .045]);
ui.Theme=Popup(traceTab,{'Dark','Light'},[.61 .145 .35 .052],'Theme',@(s,e)Dirty(h));
ui.ColorByChannel=Check(traceTab,'Color by channel (one color per channel)',false,[.04 .08 .92 .05],'ColorByChannel',@(s,e)Dirty(h));
ui.ShowLegend=Check(traceTab,'Include color legend',false,[.04 .02 .92 .05],'ShowLegend',@(s,e)Dirty(h));
set(ui.ShowLegend,'TooltipString','Molecule entries appear when their traces first enter the visible crop, then remain. Channel entries are fixed.');
set(ui.MaxDimension,'TooltipString','Canvas size before output resolution scaling. Choose 2x or 4x on the Export tab for more pixels.');
ui.IncludeBackground=Check(backgroundTab,'Include loaded background stack',~isempty(sources),[.04 .92 .92 .05], ...
    'IncludeBackground',@(s,e)Dirty(h));
names=cell(1,numel(sources));
for k=1:numel(sources)
    names{k}=sprintf('Loaded channel %d',k);if isfield(sources{k},'Name'),names{k}=sources{k}.Name;end
end
if isempty(names),names={'No stack loaded'};end
ui.BackgroundChannel=Popup(backgroundTab,names,[.04 .84 .92 .06],'BackgroundChannel',@(s,e)BackgroundChanged(h));
if isfield(m,'Channel')&&isscalar(unique([m.Channel]))
    for k=1:numel(sources)
        if isfield(sources{k},'Channel')&&sources{k}.Channel==m(1).Channel,set(ui.BackgroundChannel,'Value',k);break;end
    end
end
Label(backgroundTab,'First / last stack frame',[.04 .76 .57 .045]);
ui.BackgroundFirst=Edit(backgroundTab,'1',[.62 .76 .15 .052],'BackgroundFirst',@(s,e)Dirty(h));
ui.BackgroundLast=Edit(backgroundTab,'1',[.79 .76 .17 .052],'BackgroundLast',@(s,e)Dirty(h));
Label(backgroundTab,'Display minimum / maximum',[.04 .68 .57 .045]);
ui.BackgroundMin=Edit(backgroundTab,'0',[.62 .68 .15 .052],'BackgroundMin',@(s,e)Dirty(h));
ui.BackgroundMax=Edit(backgroundTab,'1',[.79 .68 .17 .052],'BackgroundMax',@(s,e)Dirty(h));
ui.MatchStackTime=Check(backgroundTab,'Limit movie to selected stack frames',true,[.04 .60 .92 .05],'MatchStackTime',@(s,e)Dirty(h));
ui.BackgroundInfo=Label(backgroundTab,'',[.04 .48 .92 .10]);
ui.DrawCrop=Button(backgroundTab,'Draw rectangle on aligned projection',[.04 .40 .92 .06],'DrawCrop',@(s,e)DrawCrop(h));
Label(backgroundTab,'Normalized crop [left top width height]',[.04 .33 .92 .045]);
ui.Crop=Edit(backgroundTab,'',[.04 .27 .72 .055],'Crop',@(s,e)Dirty(h));
Button(backgroundTab,'Clear',[.79 .27 .17 .055],'ClearCrop',@(s,e)ClearCrop(h));
ui.CropTiming=Popup(backgroundTab,{'Visible crop time only','Full selected trace timeline (matched inset)'}, ...
    [.04 .18 .92 .06],'CropTiming',@(s,e)Dirty(h));
Label(backgroundTab,sprintf('Draw after choosing rotation. Full timing preserves parent frames and speedups.\nThe clock is always referenced to the full selection.'),[.04 .035 .92 .12]);
Label(exportTab,'Output folder',[.04 .92 .92 .045]);
ui.Folder=Edit(exportTab,folder,[.04 .85 .72 .06],'Folder',@(s,e)Dirty(h));
Button(exportTab,'Choose...',[.79 .85 .17 .06],'ChooseFolder',@(s,e)ChooseFolder(h));
Label(exportTab,'Filename (without extension)',[.04 .77 .92 .045]);
ui.Filename=Edit(exportTab,'selected_molecule_movie',[.04 .70 .92 .06],'Filename',@(s,e)Dirty(h));
Label(exportTab,'File format',[.04 .62 .55 .045]);
ui.Format=Popup(exportTab,{'MP4','AVI','GIF'},[.61 .61 .35 .06],'Format',@(s,e)Dirty(h));
Label(exportTab,'Output resolution',[.04 .53 .55 .045]);
ui.ResolutionScale=Popup(exportTab,{'1x (standard)','2x (sharper)','4x (highest)'},[.61 .52 .35 .06],'ResolutionScale',@(s,e)Dirty(h));
set(ui.ResolutionScale,'TooltipString','Doubles or quadruples each pixel dimension, including text and trace styling. Higher settings take longer and create larger files.');
Label(exportTab,'Clock',[.04 .44 .55 .045]);
ui.TimeOrigin=Popup(exportTab,{'Elapsed from selection start','Absolute source time'},[.04 .37 .92 .06],'TimeOrigin',@(s,e)Dirty(h));
Label(exportTab,'Optional source-time range [start end] s',[.04 .28 .92 .045]);
ui.TimeRange=Edit(exportTab,'',[.04 .21 .92 .06],'TimeRange',@(s,e)Dirty(h));
Label(exportTab,sprintf('Blank time range uses the selection, stack limit and crop.\nResolution scales pixels, text and traces together.'),[.04 .045 .92 .12]);
ui.Tabs=uitabgroup(h,'Units','normalized','Position',[.39 .155 .59 .78]);
ui.ProjectionTab=uitab(ui.Tabs,'Title','Aligned projection & crop','BackgroundColor',c);
ui.MovieTab=uitab(ui.Tabs,'Title','Movie preview','BackgroundColor',c);
ui.ProjectionAxes=axes(ui.ProjectionTab,'Units','normalized','Position',[.08 .10 .88 .84], ...
    'Color',[.045 .055 .075],'XColor',fg,'YColor',fg,'YDir','reverse','DataAspectRatio',[1 1 1]);
ui.MovieAxes=axes(ui.MovieTab,'Units','normalized','Position',[.025 .14 .95 .83],'Visible','off');
ui.Slider=uicontrol(ui.MovieTab,'Style','slider','Units','normalized','Position',[.20 .055 .77 .05], ...
    'Min',1,'Max',2,'Value',1,'Tag','TraceMovieFrame','Callback',@(s,e)Preview(h));
ui.Play=Button(ui.MovieTab,'Play',[.025 .04 .15 .075],'Play',@(s,e)Play(h));
ui.Status=Label(h,'Preparing aligned projection...',[.02 .065 .63 .075]);
ui.Update=Button(h,'Update preview',[.68 .055 .14 .075],'UpdatePreview',@(s,e)Refresh(h));
ui.Export=Button(h,'Export movie',[.84 .055 .14 .075],'Export',@(s,e)Export(h));
ui.CancelButton=Button(h,'Cancel export',[.84 .01 .14 .035],'Cancel',@(s,e)Cancel(h));
set(ui.CancelButton,'Enable','off');
controls=findall(h,'Type','uicontrol');
for k=1:numel(controls)
    style=get(controls(k),'Style');set(controls(k),'ForegroundColor',fg,'FontSize',10);
    if ismember(style,{'edit','popupmenu'})
        set(controls(k),'BackgroundColor',editColor);
    elseif strcmp(style,'pushbutton')
        set(controls(k),'BackgroundColor',button);
    else
        set(controls(k),'BackgroundColor',c);
    end
end
setappdata(h,'UI',ui);BackgroundChanged(h);set(h,'Visible','on');Refresh(h);

function [m,sources,folder]=CaptureContext
global Molecule
m=Molecule;if ~isempty(m),m=m([m.Selected]==1);end
[sources,folder]=LoadedBackgrounds;

function [sources,folder]=LoadedBackgrounds
global Stack TimeInfo Config
sources={};folder=pwd;main=getappdata(0,'hMainGui');
if isempty(main),return;end
if isfield(main,'Directory')&&~isempty(main.Directory.Stack),folder=main.Directory.Stack{1};end
for k=1:numel(Stack)
    if isempty(Stack{k}),continue;end
    n=size(Stack{k},3);times=[];
    if numel(TimeInfo)>=k&&numel(TimeInfo{k})==n,times=double(TimeInfo{k}(:))/1000;end
    if isempty(times)||any(diff(times)<=0)
        dt=abs(Config.Time(min(k,numel(Config.Time))))/1000;
        if ~isfinite(dt)||dt<=0,continue;end
        times=(0:n-1)'*dt;
    end
    limits=[main.Values.ScaleMin(k) main.Values.ScaleMax(k)];
    if limits(2)<=limits(1),limits(2)=limits(1)+1;end
    b=struct('Stack',Stack{k},'TimesS',times,'PixelSizeNm',main.Values.PixSize,'Limits',limits, ...
        'Drift',[],'Name',sprintf('Loaded channel %d',k),'Channel',k);
    if n==1,b.FrameDurationS=abs(Config.Time(min(k,numel(Config.Time))))/1000;end
    if strcmp(get(main.Menu.mApplyCorrections,'Checked'),'on')
        drift=getappdata(main.fig,'Drift');
        if numel(drift)>=k&&~isempty(drift{k}),b.Drift=drift{k}(:,:,min(1:n,size(drift{k},3)));end
    end
    sources{end+1}=b; %#ok<AGROW>
end

function h=Label(parent,text,position)
h=uicontrol(parent,'Style','text','Units','normalized','Position',position,'String',text,'HorizontalAlignment','left');
function h=Edit(parent,text,position,tag,callback)
h=uicontrol(parent,'Style','edit','Units','normalized','Position',position,'String',text, ...
    'Tag',['TraceMovie' tag],'HorizontalAlignment','left','Callback',callback);
function h=Popup(parent,text,position,tag,callback)
h=uicontrol(parent,'Style','popupmenu','Units','normalized','Position',position,'String',text,'Tag',['TraceMovie' tag],'Callback',callback);
function h=Check(parent,text,value,position,tag,callback)
h=uicontrol(parent,'Style','checkbox','Units','normalized','Position',position,'String',text,'Value',value,'Tag',['TraceMovie' tag],'Callback',callback);
function h=Button(parent,text,position,tag,callback)
h=uicontrol(parent,'Style','pushbutton','Units','normalized','Position',position,'String',text,'Tag',['TraceMovie' tag],'Callback',callback);

function Dirty(h)
ui=getappdata(h,'UI');if isempty(ui)||ui.Busy,return;end
ui.Dirty=true;ui.Playing=false;set(ui.Play,'String','Play');setappdata(h,'UI',ui);
set(ui.Status,'String','Options changed. Update preview before drawing a crop or exporting.');

function RotationChanged(h)
ClearCrop(h);Refresh(h);

function BackgroundChanged(h)
ui=getappdata(h,'UI');
if isempty(ui.Sources)
    set(ui.IncludeBackground,'Value',0,'Enable','off');
    set(ui.BackgroundInfo,'String','No loaded stack. Traces can be exported on their own.');
else
    b=ui.Sources{get(ui.BackgroundChannel,'Value')};
    set(ui.BackgroundFirst,'String','1');set(ui.BackgroundLast,'String',num2str(size(b.Stack,3)));
    set(ui.BackgroundMin,'String',num2str(b.Limits(1)));set(ui.BackgroundMax,'String',num2str(b.Limits(2)));
    set(ui.BackgroundInfo,'String',sprintf('%g nm/pixel; %.3f–%.3f s. Stack follows source time; outside its range the background is blank.', ...
        b.PixelSizeNm,b.TimesS(1),b.TimesS(end)));
end
Dirty(h);

function values=Numbers(control)
text=strtrim(get(control,'String'));
if isempty(text),values=[];return;end
text=regexprep(text,'[\[\],;]',' ');parts=regexp(strtrim(text),'\s+','split');values=str2double(parts);
if any(~isfinite(values)),error('FIESTA:TraceMovieInput','Enter finite numbers in %s.',get(control,'Tag'));end

function o=ReadOptions(h)
ui=getappdata(h,'UI');o=fTraceMovie('defaults');
keys={'PlaybackSpeed','FrameRate','DecayMs','BaselineAlpha','MarkerSize','LineWidth','DeadTimeSpeedup','BufferFraction','MaxDimension','ScaleBarNm'};
for k=1:numel(keys),o.(keys{k})=Numbers(ui.(keys{k}));end
rotations={'none','horizontal','vertical'};themes={'dark','light'};origins={'relative','absolute'};timing={'visible','full'};
o.Rotation=rotations{get(ui.Rotation,'Value')};o.Theme=themes{get(ui.Theme,'Value')};
o.TimeOrigin=origins{get(ui.TimeOrigin,'Value')};o.ShowLegend=logical(get(ui.ShowLegend,'Value'));
if isfield(ui,'ColorByChannel'),o.ColorByChannel=logical(get(ui.ColorByChannel,'Value'));end
if isfield(ui,'ResolutionScale'),scales=[1 2 4];o.ResolutionScale=scales(get(ui.ResolutionScale,'Value'));end
o.CropNormalized=Numbers(ui.Crop);o.CropTiming=timing{get(ui.CropTiming,'Value')};o.TimeRangeS=Numbers(ui.TimeRange);
if get(ui.IncludeBackground,'Value')
    b=ui.Sources{get(ui.BackgroundChannel,'Value')};n=size(b.Stack,3);
    first=Numbers(ui.BackgroundFirst);last=Numbers(ui.BackgroundLast);
    validateattributes(first,{'numeric'},{'scalar','integer','>=',1,'<=',n});
    validateattributes(last,{'numeric'},{'scalar','integer','>=',first,'<=',n});
    b.Stack=b.Stack(:,:,first:last);b.TimesS=b.TimesS(first:last);
    if isfield(b,'Drift')&&size(b.Drift,3)>1,b.Drift=b.Drift(:,:,first:last);end
    b.Limits=[Numbers(ui.BackgroundMin) Numbers(ui.BackgroundMax)];
    if ~isfield(b,'FrameDurationS')||isempty(b.FrameDurationS)
        source=ui.Sources{get(ui.BackgroundChannel,'Value')};b.FrameDurationS=1;
        if numel(source.TimesS)>1,b.FrameDurationS=median(diff(source.TimesS));end
    end
    o.Background=b;
    if get(ui.MatchStackTime,'Value')
        limits=b.TimesS([1 end])';
        if isempty(o.TimeRangeS),o.TimeRangeS=limits;
        else
            validateattributes(o.TimeRangeS,{'numeric'},{'vector','numel',2});
            o.TimeRangeS=[max(limits(1),o.TimeRangeS(1)) min(limits(2),o.TimeRangeS(2))];
        end
    end
end

function Refresh(h)
ui=getappdata(h,'UI');if ui.Busy,return;end
ui.Playing=false;set(ui.Play,'String','Play');set(ui.Status,'String','Updating aligned projection and movie preview...');drawnow;
try
    o=ReadOptions(h);baseOptions=o;baseOptions.CropNormalized=[];
    base=fTraceMovie('prepare',ui.Molecules,baseOptions);
    if isempty(o.CropNormalized),plan=base;else,plan=fTraceMovie('prepare',ui.Molecules,o);end
    if ~isempty(ui.Renderer),fTraceMovie('dispose',ui.Renderer);end
    ui.Renderer=fTraceMovie('renderer',plan);ui.Plan=plan;ui.BasePlan=base;ui.Dirty=false;
    ui.ROI=[];setappdata(h,'UI',ui);cla(ui.ProjectionAxes);hold(ui.ProjectionAxes,'on');
    if ~isempty(o.Background)
        sz=round(base.PlotRectPx([4 3]));sz=max(1,round(sz*min(1,600/max(sz))));
        bg=fTraceMovieBackground('prepare',o.Background,base,sz);gray=fTraceMovieBackground('projection',bg,true);
        span=diff(base.BoundsNm,1,1);dx=span(1)/sz(2);dy=span(2)/sz(1);
        image(ui.ProjectionAxes,'XData',base.BoundsNm(:,1)'+[dx/2 -dx/2], ...
            'YData',base.BoundsNm(:,2)'+[dy/2 -dy/2],'CData',repmat(gray,1,1,3));
    end
    for j=1:base.MoleculeCount
        xy=base.Traces(j).XYNm;plot(ui.ProjectionAxes,xy(:,1),xy(:,2),'-','Color',base.Colors(j,:),'LineWidth',.5);
    end
    set(ui.ProjectionAxes,'XLim',base.BoundsNm(:,1)','YLim',base.BoundsNm(:,2)','YDir','reverse','DataAspectRatio',[1 1 1]);
    xlabel(ui.ProjectionAxes,'Aligned X (nm)');ylabel(ui.ProjectionAxes,'Aligned Y (nm)');
    if ~isempty(o.CropNormalized)
        bounds=plan.BoundsNm;ui.ROI=drawrectangle(ui.ProjectionAxes,'Position',[bounds(1,:) diff(bounds,1,1)], ...
            'Color',[.2 .9 1],'FixedAspectRatio',false,'Selected',false, ...
            'DrawingArea',[base.BoundsNm(1,:) diff(base.BoundsNm,1,1)]);
        addlistener(ui.ROI,'ROIMoved',@(s,e)StoreCrop(h));
    end
    maxFrame=max(2,plan.FrameCount);value=min(maxFrame,max(1,round(get(ui.Slider,'Value'))));
    set(ui.Slider,'Min',1,'Max',maxFrame,'Value',value,'SliderStep',min(1,[1 10]/max(1,plan.FrameCount-1)));
    setappdata(h,'UI',ui);Preview(h);
    set(ui.Status,'String',sprintf('%d molecules; %d x %d px; %d output images; %.2f s playback; source %.3f–%.3f s. Legend %s.', ...
        plan.MoleculeCount,plan.FrameSize,plan.FrameCount,plan.FrameCount/o.FrameRate,plan.TimeRangeS,OnOff(o.ShowLegend)));
catch err
    setappdata(h,'LastError',err);
    set(ui.Status,'String',err.message);ui=getappdata(h,'UI');ui.Dirty=true;setappdata(h,'UI',ui);
end

function text=OnOff(value)
text='off';if value,text='on';end

function Preview(h)
ui=getappdata(h,'UI');if isempty(ui.Plan)||ui.Dirty,return;end
k=min(ui.Plan.FrameCount,max(1,round(get(ui.Slider,'Value'))));
rgb=fTraceMovie('render',ui.Renderer,ui.Plan,k);image(ui.MovieAxes,rgb);axis(ui.MovieAxes,'image','off');

function Play(h)
ui=getappdata(h,'UI');
if ui.Playing,ui.Playing=false;setappdata(h,'UI',ui);return;end
if ui.Dirty,Refresh(h);ui=getappdata(h,'UI');end
if ui.Dirty,return;end
ui.Playing=true;set(ui.Play,'String','Pause');setappdata(h,'UI',ui);
ui.Tabs.SelectedTab=ui.MovieTab;start=round(get(ui.Slider,'Value'));if start>=ui.Plan.FrameCount,start=1;end
tic;previous=0;
while isgraphics(h)
    ui=getappdata(h,'UI');if ~ui.Playing||ui.Dirty||ui.Busy,break;end
    k=min(ui.Plan.FrameCount,start+floor(toc*ui.Plan.Options.FrameRate));
    if k~=previous,set(ui.Slider,'Value',k);Preview(h);previous=k;else,drawnow;pause(.005);end
    if k==ui.Plan.FrameCount,break;end
end
if isgraphics(h),ui=getappdata(h,'UI');ui.Playing=false;set(ui.Play,'String','Play');setappdata(h,'UI',ui);end

function DrawCrop(h)
ui=getappdata(h,'UI');if ui.Busy,return;end
if ui.Dirty,Refresh(h);ui=getappdata(h,'UI');end
if ui.Dirty,return;end
ui.Tabs.SelectedTab=ui.ProjectionTab;
if ~isempty(ui.ROI)&&isvalid(ui.ROI),delete(ui.ROI);end
set(ui.Status,'String','Drag a rectangular crop on the aligned projection. Move or resize either dimension afterward.');
ui.ROI=drawrectangle(ui.ProjectionAxes,'Color',[.2 .9 1],'FixedAspectRatio',false, ...
    'DrawingArea',[ui.BasePlan.BoundsNm(1,:) diff(ui.BasePlan.BoundsNm,1,1)]);
setappdata(h,'UI',ui);addlistener(ui.ROI,'ROIMoved',@(s,e)StoreCrop(h));StoreCrop(h);

function StoreCrop(h)
ui=getappdata(h,'UI');if isempty(ui.ROI)||~isvalid(ui.ROI),return;end
position=ui.ROI.Position;
if isempty(position)||any(position(3:4)<=0),return;end
bounds=ui.BasePlan.FullBoundsNm;span=diff(bounds,1,1);
box=[(position(1:2)-bounds(1,:))./span position(3:4)./span];
box=max(0,min(1,box));box(3:4)=min(box(3:4),1-box(1:2));
set(ui.Crop,'String',sprintf('%.10g ',box));Dirty(h);
set(ui.Status,'String','Crop stored. Update preview to apply its framing and time range.');

function ClearCrop(h)
ui=getappdata(h,'UI');if ~isempty(ui.ROI)&&isvalid(ui.ROI),delete(ui.ROI);end
ui.ROI=[];set(ui.Crop,'String','');setappdata(h,'UI',ui);Dirty(h);

function ChooseFolder(h)
ui=getappdata(h,'UI');folder=uigetdir(get(ui.Folder,'String'),'Choose movie output folder');
if isequal(folder,0),return;end
set(ui.Folder,'String',folder);setappdata(0,'TraceMovieOutputFolder',folder);

function Export(h)
ui=getappdata(h,'UI');if ui.Busy,return;end
if ui.Dirty,Refresh(h);ui=getappdata(h,'UI');end
if ui.Dirty,return;end
try
    o=ReadOptions(h);folder=strtrim(get(ui.Folder,'String'));name=strtrim(get(ui.Filename,'String'));
    [pathPart,~,extension]=fileparts(name);
    if isempty(name)||~isempty(pathPart)||~isempty(extension),error('FIESTA:TraceMovieName','Enter a filename without a folder or extension.');end
    if ~isfolder(folder),error('FIESTA:TraceMovieFolder','Choose an existing output folder.');end
    formats={'.mp4','.avi','.gif'};file=fullfile(folder,[name formats{get(ui.Format,'Value')}]);
    if isfile(file)
        if ~strcmp(questdlg(sprintf('Replace %s?',file),'Movie exists','Replace','Cancel','Cancel'),'Replace'),return;end
        o.Overwrite=true;
    end
    ui.Busy=true;ui.Cancel=false;ui.Playing=false;setappdata(h,'UI',ui);Busy(h,true);
    if ~isempty(ui.ROI)&&isvalid(ui.ROI),ui.ROI.Selected=false;ui.ROI.InteractionsAllowed='none';drawnow;end
    o.ProgressFcn=@(done,total)Progress(h,done,total);
    report=fTraceMovie('export',ui.Molecules,file,o);
    if isgraphics(h)
        setappdata(h,'LastReport',report);setappdata(0,'TraceMovieOutputFolder',folder);
        set(ui.Status,'String',sprintf('Saved %s (%.2f s playback).',file,report.MovieDurationS));
    end
catch err
    if isgraphics(h),set(ui.Status,'String',err.message);end
end
if isgraphics(h)
    ui=getappdata(h,'UI');ui.Busy=false;setappdata(h,'UI',ui);Busy(h,false);
    if ~isempty(ui.ROI)&&isvalid(ui.ROI),ui.ROI.InteractionsAllowed='all';end
    if ui.Closing,Close(h);end
end

function Busy(h,value)
ui=getappdata(h,'UI');controls=findall(h,'Type','uicontrol');set(controls,'Enable',OnOff(~value));
set(ui.CancelButton,'Enable',OnOff(value));
if isempty(ui.Sources),set(ui.IncludeBackground,'Enable','off');end

function proceed=Progress(h,done,total)
proceed=false;if ~isgraphics(h),return;end
ui=getappdata(h,'UI');
set(ui.Status,'String',sprintf('Writing output image %d / %d (%.0f%%).',done,total,100*done/total));drawnow;
if isgraphics(h),ui=getappdata(h,'UI');proceed=~ui.Cancel;end

function Cancel(h)
ui=getappdata(h,'UI');ui.Cancel=true;setappdata(h,'UI',ui);

function Close(h)
ui=getappdata(h,'UI');
if ~isempty(ui)&&ui.Busy,ui.Cancel=true;ui.Closing=true;setappdata(h,'UI',ui);return;end
if ~isempty(ui)&&~isempty(ui.Renderer),fTraceMovie('dispose',ui.Renderer);end
delete(h);
