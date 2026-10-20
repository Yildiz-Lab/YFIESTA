function varargout=fSpatialBackgroundGui(action,varargin)
%FSPATIALBACKGROUNDGUI Editor, threshold preview and advisory warning delivery.
% Actions: create/commit/preview/warnings. Configuration OK commits edits.
global Config Stack
switch action
    case 'create'
        parent=varargin{1};
        setappdata(parent.fig,'SpatialBackgroundSettings',Config.Threshold);
        set(parent.pGeneral.rFilterNone,'Position',[.025 .1 .12 .8]);
        set(parent.pGeneral.rFilterAverage,'Position',[.15 .1 .20 .8]);
        set(parent.pGeneral.rFilterSmooth,'Position',[.36 .1 .18 .8]);
        set(parent.pGeneral.cFilterBackground,'Position',[.55 .1 .20 .8]);
        button=uicontrol(parent.pGeneral.pFilter,'Style','pushbutton','Units','normalized', ...
            'Position',[.77 .15 .21 .65],'String','Spatial background...', ...
            'BackgroundColor',[.16 .16 .16],'ForegroundColor',[.9 .9 .9], ...
            'Tag','SpatialBackgroundSettings','FontSize',10, ...
            'Callback',@(src,event) OpenEditor(parent));
        setappdata(parent.fig,'SpatialBackgroundButton',button);
    case 'commit'
        parent=varargin{1}; threshold=varargin{2};
        edited=getappdata(parent.fig,'SpatialBackgroundSettings');
        if isfield(edited,'SpatialBackground')
            threshold.SpatialBackground=edited.SpatialBackground;
        end
        for channel=1:max(1,numel(Stack))
            s=fSpatialBackground('get',threshold,channel);
            imageSize=[];
            if channel<=numel(Stack), imageSize=[size(Stack{channel},1) size(Stack{channel},2)]; end
            fSpatialBackground('validate',s,imageSize,threshold.Mode);
            if s.Enabled && contains(threshold.Filter,'+background')
                error('FIESTA:BackgroundConflict','Turn off filter background when using a spatial background image.');
            end
        end
        varargout{1}=threshold;
    case 'preview'
        varargout{1}=ThresholdPreview(varargin{:});
    case 'warnings'
        varargout{1}=CheckWarnings(varargin{:});
end

function OpenEditor(parent)
global Stack
existing=findall(0,'Tag','SpatialBackgroundDialog');
for k=1:numel(existing)
    previous=getappdata(existing(k),'UI');
    if ~isempty(previous) && isequal(previous.Parent.fig,parent.fig)
        figure(existing(k)); return;
    end
end
main=getappdata(0,'hMainGui'); count=max(1,numel(Stack));
c=get(parent.fig,'Color'); ctext=.9*ones(1,3); ctextbox=.12*ones(1,3); cbutton=.16*ones(1,3);
h=dialog('Name','Spatial background for thresholding','Tag','SpatialBackgroundDialog', ...
    'Units','pixels','Position',[100 100 1100 660],'Resize','on','WindowStyle','normal', ...
    'Color',c,'DefaultUicontrolBackgroundColor',c,'DefaultUicontrolForegroundColor',ctext);
movegui(h,'center');
ui.Parent=parent; ui.Threshold=getappdata(parent.fig,'SpatialBackgroundSettings');
ui.Channel=min(count,main.Values.FrameIdx(1)); ui.Drafts=cell(1,count);
left=uipanel(h,'Units','normalized','Position',[.02 .17 .36 .78],'BorderType','none','BackgroundColor',c);
uicontrol(left,'Style','text','Units','normalized','Position',[.02 .93 .24 .045], ...
    'String','Channel','HorizontalAlignment','left');
ui.ChannelBox=uicontrol(left,'Style','popupmenu','Units','normalized','Position',[.28 .92 .70 .065], ...
    'String',arrayfun(@(n)sprintf('Channel %d',n),1:count,'UniformOutput',false), ...
    'Value',ui.Channel,'Tag','SpatialBackgroundChannel','Callback',@(src,event) ChangeChannel(h));
source=uipanel(left,'Units','normalized','Position',[0 .48 1 .40],'Title','Source', ...
    'BackgroundColor',c,'ForegroundColor',ctext);
ui.Source=uicontrol(source,'Style','text','Units','normalized','Position',[.04 .66 .92 .26], ...
    'HorizontalAlignment','left','Tag','SpatialBackgroundSource');
uicontrol(source,'Style','pushbutton','Units','normalized','Position',[.04 .48 .62 .14], ...
    'String','Choose TIFF...','Tag','SpatialBackgroundLoad','Callback',@(src,event) LoadMap(h));
uicontrol(source,'Style','pushbutton','Units','normalized','Position',[.70 .48 .26 .14], ...
    'String','Clear','Tag','SpatialBackgroundClear','Callback',@(src,event) ClearMap(h));
ui.Direct=uicontrol(source,'Style','radiobutton','Units','normalized','Position',[.04 .32 .45 .13], ...
    'String','Use TIFF directly','Tag','SpatialBackgroundDirect','Callback',@(src,event) SourceMode(h,false));
ui.Gaussian=uicontrol(source,'Style','radiobutton','Units','normalized','Position',[.51 .32 .45 .13], ...
    'String','Gaussian smooth','Tag','SpatialBackgroundGaussian','Callback',@(src,event) SourceMode(h,true));
uicontrol(source,'Style','text','Units','normalized','Position',[.04 .18 .61 .10], ...
    'String','Gaussian sigma (pixels)','HorizontalAlignment','left');
ui.Sigma=uicontrol(source,'Style','edit','Units','normalized','Position',[.70 .17 .26 .14], ...
    'String','30','Tag','SpatialBackgroundSigma','Callback',@(src,event) ControlsChanged(h), ...
    'TooltipString','Positive Gaussian standard deviation in pixels. Generate again after changing sigma.');
ui.CreateBlur=uicontrol(source,'Style','pushbutton','Units','normalized','Position',[.04 .02 .92 .14], ...
    'String','Generate & save background','Tag','SpatialBackgroundCreateBlur','Callback',@(src,event) CreateBlur(h));
correction=uipanel(left,'Units','normalized','Position',[0 .07 1 .38],'Title','Threshold correction', ...
    'BackgroundColor',c,'ForegroundColor',ctext);
ui.Enabled=uicontrol(correction,'Style','checkbox','Units','normalized','Position',[.04 .80 .92 .13], ...
    'String','Enable spatial background correction','Tag','SpatialBackgroundEnabled','Callback',@(src,event) ControlsChanged(h));
ui.File=uicontrol(correction,'Style','text','Units','normalized','Position',[.04 .58 .92 .20], ...
    'HorizontalAlignment','left','Tag','SpatialBackgroundFile','FontSize',9);
uicontrol(correction,'Style','text','Units','normalized','Position',[.04 .44 .61 .10], ...
    'String','Camera baseline (counts)','HorizontalAlignment','left');
ui.Baseline=uicontrol(correction,'Style','edit','Units','normalized','Position',[.70 .43 .26 .14], ...
    'Tag','SpatialBackgroundBaseline','Callback',@(src,event) ControlsChanged(h));
uicontrol(correction,'Style','text','Units','normalized','Position',[.04 .27 .61 .10], ...
    'String','Correction strength (1 = full)','HorizontalAlignment','left');
ui.Strength=uicontrol(correction,'Style','edit','Units','normalized','Position',[.70 .26 .26 .14], ...
    'Tag','SpatialBackgroundStrength','Callback',@(src,event) ControlsChanged(h));
ui.Preview=uicontrol(correction,'Style','pushbutton','Units','normalized','Position',[.04 .04 .92 .15], ...
    'String','Preview threshold mask','Tag','SpatialBackgroundPreview','Callback',@(src,event) Preview(h));
uicontrol(left,'Style','text','Units','normalized','Position',[.02 0 .96 .065],'FontSize',9, ...
    'String','1x scales spatial variation to each frame above the baseline.','HorizontalAlignment','left');
ui.TabGroup=uitabgroup(h,'Units','normalized','Position',[.40 .19 .58 .76], ...
    'Tag','SpatialBackgroundTabs');
keys={'Reference','Background','Threshold'}; labels={'Reference','Background','Threshold mask'};
for k=1:numel(keys)
    key=keys{k};
    ui.Tabs.(key)=uitab(ui.TabGroup,'Title',labels{k},'BackgroundColor',c,'Tag',['SpatialBackgroundTab' key]);
    ui.Axes.(key)=axes('Parent',ui.Tabs.(key),'Units','normalized','Position',[.12 .19 .72 .72], ...
        'Color',ctextbox,'XColor',ctext,'YColor',ctext,'Visible','off','Tag',['SpatialBackgroundAxes' key]);
    ui.Empty.(key)=uicontrol(ui.Tabs.(key),'Style','text','Units','normalized', ...
        'Position',[.06 .38 .88 .20],'HorizontalAlignment','center','FontSize',11);
    ui.Info.(key)=uicontrol(ui.Tabs.(key),'Style','text','Units','normalized', ...
        'Position',[.04 .015 .92 .11],'HorizontalAlignment','left','FontSize',9);
end
ui.Status=uicontrol(h,'Style','text','Units','normalized','Position',[.02 .09 .36 .065], ...
    'HorizontalAlignment','left','Tag','SpatialBackgroundStatus','FontSize',9);
ui.Warning=uicontrol(h,'Style','text','Units','normalized','Position',[.41 .09 .56 .075], ...
    'HorizontalAlignment','left','ForegroundColor',[1 .7 .3],'Tag','SpatialBackgroundWarning','FontSize',9);
uicontrol(h,'Style','text','Units','normalized','Position',[.02 .025 .46 .045], ...
    'String','Changes take effect when Configuration is applied.','HorizontalAlignment','left','FontSize',9);
uicontrol(h,'Style','pushbutton','Units','normalized','Position',[.51 .025 .30 .05], ...
    'String','Apply to Configuration','Tag','SpatialBackgroundApply','Callback',@(src,event) Apply(h));
uicontrol(h,'Style','pushbutton','Units','normalized','Position',[.83 .025 .15 .05], ...
    'String','Cancel','Callback',@(src,event) delete(h));
set(findall(h,'Style','pushbutton'),'BackgroundColor',cbutton);
set([ui.ChannelBox ui.Baseline ui.Strength ui.Sigma],'BackgroundColor',ctextbox);
setappdata(h,'UI',ui); Populate(h);
set(ui.TabGroup,'SelectionChangedFcn',@(src,event) TabSelected(h));

function d=InitialDraft(s)
source=s.SourceFile; if isempty(source), source=s.File; end
sigma='30'; if ~isempty(s.BlurSigma), sigma=num2str(s.BlurSigma,15); end
d=struct('Source',source,'Gaussian',~isempty(s.BlurSigma),'Sigma',sigma, ...
    'Reference',[],'Background',[],'Normalized',false);
% Stored maps remain usable even when the original TIFF has moved.
try
    if isfile(source), d.Reference=fSpatialBackground('reference',source); end
    if strcmp(source,s.File) && ~isempty(d.Reference)
        d.Background=d.Reference;
    elseif isfile(s.File)
        d.Background=fSpatialBackground('reference',s.File);
    end
catch
    % The missing reference is reported inline; the saved map is portable.
end
if isempty(d.Background), d.Background=s.NormalizedMap; d.Normalized=true; end

function Populate(h)
u=getappdata(h,'UI'); s=fSpatialBackground('get',u.Threshold,u.Channel);
if isempty(u.Drafts{u.Channel}), u.Drafts{u.Channel}=InitialDraft(s); end
d=u.Drafts{u.Channel}; setappdata(h,'UI',u);
set(u.Enabled,'Value',s.Enabled); set(u.Baseline,'String',num2str(s.Baseline)); set(u.Strength,'String',num2str(s.Strength));
set(u.Direct,'Value',~d.Gaussian); set(u.Gaussian,'Value',d.Gaussian); set(u.Sigma,'String',d.Sigma);
enable='off'; if d.Gaussian, enable='on'; end
set(u.Sigma,'Enable',enable);
if isempty(d.Source), enable='off'; end
set(u.CreateBlur,'Enable',enable);
if isempty(d.Source)
    sourceText='No reference selected.';
else
    sourceText=FileLabel(d.Source);
end
set(u.Source,'String',sourceText,'TooltipString',d.Source);
RenderImage(h,'Reference',d.Reference,[],'Reference (mean of TIFF pages)');
set(u.Info.Reference,'String','TIFF pages are averaged. Display contrast does not change saved intensity counts.');
RefreshBackground(h,s);
InvalidateMask(h); UpdateWarning(h,s);
[ready,message]=BackgroundReady(u,s,u.Channel);
if ready, SetStatus(h,'Background ready. Preview the threshold mask before applying.'); else, SetStatus(h,message); end

function RefreshBackground(h,s)
u=getappdata(h,'UI'); d=u.Drafts{u.Channel}; [ready,message]=BackgroundReady(u,s,u.Channel);
if ready
    text=['Background: ' FileLabel(s.File)];
    limits=DisplayLimits(d.Reference);
    title='Background (counts)';
    if d.Normalized, limits=[]; title='Stored normalized background (median = 1)'; end
    RenderImage(h,'Background',d.Background,limits,title);
    set(u.Info.Background,'String',['Selected background: ' s.File]);
else
    text=message;
    EmptyPreview(h,'Background',message);
end
set(u.File,'String',text,'TooltipString',s.File);

function s=Collect(h)
u=getappdata(h,'UI'); s=fSpatialBackground('get',u.Threshold,u.Channel);
baseline=str2double(get(u.Baseline,'String')); strength=str2double(get(u.Strength,'String'));
validateattributes(baseline,{'numeric'},{'real','scalar','finite','nonnegative'});
validateattributes(strength,{'numeric'},{'real','scalar','finite','nonnegative'});
if baseline~=s.Baseline && ~isempty(s.NormalizedMap)
    reloaded=fSpatialBackground('load',s.File,baseline,strength);
    reloaded.SourceFile=s.SourceFile; reloaded.BlurSigma=s.BlurSigma; s=reloaded;
end
s.Baseline=baseline; s.Strength=strength; s.Enabled=logical(get(u.Enabled,'Value'));

function [u,s]=Capture(h)
u=getappdata(h,'UI'); s=Collect(h);
u.Threshold=fSpatialBackground('set',u.Threshold,u.Channel,s);
d=u.Drafts{u.Channel}; d.Gaussian=logical(get(u.Gaussian,'Value')); d.Sigma=get(u.Sigma,'String');
u.Drafts{u.Channel}=d; setappdata(h,'UI',u);

function [ready,message]=BackgroundReady(u,s,channel)
d=u.Drafts{channel}; ready=false;
message='Choose a TIFF to get started.';
if isempty(d) || isempty(d.Source)
    if ~isempty(s.NormalizedMap), ready=true; message=''; end
    return;
end
if d.Gaussian
    sigma=str2double(d.Sigma);
    if ~isfinite(sigma) || ~isreal(sigma) || sigma<=0 || isempty(s.BlurSigma) || ...
            sigma~=s.BlurSigma || ~strcmp(d.Source,s.SourceFile)
        message='Click Generate & save background to use the selected source and sigma.'; return;
    end
elseif ~strcmp(d.Source,s.File) || ~isempty(s.BlurSigma)
    message='Choose Use TIFF directly to use the selected reference as the background.'; return;
end
ready=~isempty(s.NormalizedMap);
if ready, message=''; end

function RequireReady(u,s,channel)
[ready,message]=BackgroundReady(u,s,channel);
if ~ready, error('FIESTA:BackgroundNotReady','Channel %d: %s',channel,message); end

function ChangeChannel(h)
u=getappdata(h,'UI'); previous=u.Channel;
try
    [u,~]=Capture(h); u.Channel=get(u.ChannelBox,'Value'); setappdata(h,'UI',u); Populate(h);
catch err
    set(u.ChannelBox,'Value',previous); SetStatus(h,err.message,true);
end

function SourceMode(h,gaussian)
u=getappdata(h,'UI'); set(u.Gaussian,'Value',gaussian); set(u.Direct,'Value',~gaussian);
try
    [u,s]=Capture(h); d=u.Drafts{u.Channel};
    if ~gaussian && ~isempty(d.Source)
        [s,d.Background]=fSpatialBackground('load',d.Source,s.Baseline,s.Strength);
        ValidateImageSize(s,u.Channel); d.Normalized=false;
        u.Threshold=fSpatialBackground('set',u.Threshold,u.Channel,s);
        u.Drafts{u.Channel}=d; setappdata(h,'UI',u);
    end
    Populate(h);
catch err
    InvalidateMask(h); SetStatus(h,err.message,true);
end

function ControlsChanged(h)
u=getappdata(h,'UI'); d=u.Drafts{u.Channel}; d.Sigma=get(u.Sigma,'String');
u.Drafts{u.Channel}=d; setappdata(h,'UI',u); InvalidateMask(h);
s=fSpatialBackground('get',u.Threshold,u.Channel); [ready,message]=BackgroundReady(u,s,u.Channel);
RefreshBackground(h,s);
if ~ready
    SetStatus(h,message);
else
    SetStatus(h,'Settings changed. Preview the threshold mask to refresh.');
end
set(u.Warning,'String','');

function LoadMap(h)
try
    [u,s]=Capture(h);
    [file,folder]=uigetfile({'*.tif;*.tiff','Reference TIFF (*.tif, *.tiff)'},'Choose spatial background reference');
    if isequal(file,0), return; end
    file=fullfile(folder,file); d=u.Drafts{u.Channel}; d.Source=file;
    if d.Gaussian
        d.Reference=fSpatialBackground('reference',file); d.Background=[]; d.Normalized=false;
    else
        [s,d.Reference]=fSpatialBackground('load',file,s.Baseline,s.Strength);
        ValidateImageSize(s,u.Channel); d.Background=d.Reference; d.Normalized=false;
        u.Threshold=fSpatialBackground('set',u.Threshold,u.Channel,s);
    end
    u.Drafts{u.Channel}=d; setappdata(h,'UI',u); Populate(h);
    set(u.TabGroup,'SelectedTab',u.Tabs.Reference);
catch err
    SetStatus(h,err.message,true);
end

function ValidateImageSize(s,channel)
global Stack
imageSize=[];
if channel<=numel(Stack), imageSize=[size(Stack{channel},1) size(Stack{channel},2)]; end
fSpatialBackground('validate',s,imageSize,'constant');

function CreateBlur(h)
global Stack
try
    [u,s]=Capture(h); d=u.Drafts{u.Channel};
    if isempty(d.Source), error('FIESTA:NoReference','Choose a reference TIFF first.'); end
    imageSize=[];
    if u.Channel<=numel(Stack), imageSize=[size(Stack{u.Channel},1) size(Stack{u.Channel},2)]; end
    [s,d.Reference,d.Background]=fSpatialBackground('gaussian',d.Source,str2double(d.Sigma),s.Baseline,s.Strength,imageSize);
    d.Normalized=false; d.Gaussian=true; d.Sigma=num2str(s.BlurSigma,15);
    u.Threshold=fSpatialBackground('set',u.Threshold,u.Channel,s); u.Drafts{u.Channel}=d;
    setappdata(h,'UI',u); Populate(h); set(u.TabGroup,'SelectedTab',u.Tabs.Background);
    SetStatus(h,['Saved background: ' s.File]);
catch err
    SetStatus(h,err.message,true);
end

function ClearMap(h)
u=getappdata(h,'UI'); u.Threshold=fSpatialBackground('set',u.Threshold,u.Channel,fSpatialBackground('defaults'));
u.Drafts{u.Channel}=[]; setappdata(h,'UI',u); Populate(h);

function Apply(h)
try
    [u,~]=Capture(h);
    mode='constant'; if get(u.Parent.pGeneral.rVariable,'Value'), mode='variable'; end
    for channel=1:numel(u.Drafts)
        s=fSpatialBackground('get',u.Threshold,channel);
        if ~s.Enabled, continue; end
        if ~isempty(u.Drafts{channel}), RequireReady(u,s,channel); end
        ValidateImageSize(s,channel); fSpatialBackground('validate',s,[],mode);
    end
    setappdata(u.Parent.fig,'SpatialBackgroundSettings',u.Threshold);
    if fSpatialBackground('enabled',u.Threshold,1:numel(u.Drafts))
        set(u.Parent.pGeneral.cFilterBackground,'Value',0);
    end
    delete(h);
catch err
    SetStatus(h,err.message,true);
end

function Preview(h)
global Stack
try
    [u,s]=Capture(h); RequireReady(u,s,u.Channel);
    if ~s.Enabled, error('FIESTA:SpatialBackgroundDisabled','Enable spatial background correction to preview the mask.'); end
    main=getappdata(0,'hMainGui'); channel=u.Channel;
    if channel>numel(Stack) || isempty(Stack{channel})
        error('FIESTA:NoImage','Load a movie to preview its threshold mask. Reference and Background previews remain available.');
    end
    mode='constant'; if get(u.Parent.pGeneral.rRelative,'Value'), mode='relative'; end
    if get(u.Parent.pGeneral.rVariable,'Value'), mode='variable'; end
    fSpatialBackground('validate',s,[size(Stack{channel},1) size(Stack{channel},2)],mode);
    frame=main.Values.FrameIdx(min(channel+1,numel(main.Values.FrameIdx)));
    if ~isreal(frame) || frame<1, frame=1; end
    frame=min(frame,size(Stack{channel},3)); image=double(Stack{channel}(:,:,frame));
    filter='none';
    if get(u.Parent.pGeneral.rFilterAverage,'Value'), filter='average'; end
    if get(u.Parent.pGeneral.rFilterSmooth,'Value'), filter='smooth'; end
    p=struct('spatial_background',s,'binary_image_processing',filter,'background_filter','', ...
        'threshold',main.Values.Thresh(channel));
    if strcmp(mode,'relative'), p.threshold=main.Values.RelThresh(channel)*1i; end
    [mask,details]=fSpatialBackground('threshold',image,p);
    RenderImage(h,'Threshold',mask,[0 1],sprintf('Detection mask - channel %d, frame %d',channel,frame));
    set(u.Info.Threshold,'String',sprintf('Native image pixels | cutoff %.6g | strength %.6g | %s filter', ...
        details.Threshold,s.Strength,filter));
    setappdata(h,'ThresholdPreviewDetails',details); setappdata(h,'ThresholdPreviewMask',mask);
    UpdateWarning(h,s,image); SetStatus(h,'Threshold preview is current. Raw movie pixels are unchanged.');
    % Temporarily suppress the selection callback to avoid calculating twice.
    callback=get(u.TabGroup,'SelectionChangedFcn'); set(u.TabGroup,'SelectionChangedFcn',[]);
    set(u.TabGroup,'SelectedTab',u.Tabs.Threshold); set(u.TabGroup,'SelectionChangedFcn',callback);
catch err
    InvalidateMask(h); SetStatus(h,err.message,true);
end

function TabSelected(h)
u=getappdata(h,'UI');
if isequal(get(u.TabGroup,'SelectedTab'),u.Tabs.Threshold), Preview(h); end

function limits=DisplayLimits(image)
limits=[]; if isempty(image), return; end
values=sort(double(image(:))); count=numel(values);
limits=[values(max(1,round(.01*count))) values(max(1,round(.99*count)))];
if limits(1)==limits(2), limits=limits+[-.5 .5]; end

function text=FileLabel(file)
[~,name,ext]=fileparts(file); text=[name ext];
if numel(text)>80, text=[text(1:36) ' ... ' text(end-38:end)]; end

function RenderImage(h,key,image,limits,label)
u=getappdata(h,'UI');
if isempty(image), EmptyPreview(h,key,'No image available. Choose a reference TIFF.'); return; end
if isempty(limits), limits=DisplayLimits(image); end
ax=u.Axes.(key); imagesc(ax,image,limits); axis(ax,'image');
set(ax,'Visible','on','Color',[.12 .12 .12],'XColor',[.9 .9 .9],'YColor',[.9 .9 .9]);
title(ax,label,'Color',[.9 .9 .9],'Interpreter','none'); colormap(ax,'gray');
bar=colorbar(ax); set(bar,'Color',[.9 .9 .9]); set(u.Empty.(key),'Visible','off');

function EmptyPreview(h,key,message)
u=getappdata(h,'UI'); ax=u.Axes.(key); cla(ax); title(ax,''); colorbar(ax,'off'); set(ax,'Visible','off');
set(u.Empty.(key),'String',message,'Visible','on'); set(u.Info.(key),'String','');

function InvalidateMask(h)
EmptyPreview(h,'Threshold','Click Preview threshold mask to show the current movie frame.');
if isappdata(h,'ThresholdPreviewMask'), rmappdata(h,'ThresholdPreviewMask'); end
if isappdata(h,'ThresholdPreviewDetails'), rmappdata(h,'ThresholdPreviewDetails'); end

function SetStatus(h,message,isError)
u=getappdata(h,'UI'); color=[.9 .9 .9];
if nargin>2 && isError, color=[1 .55 .4]; end
set(u.Status,'String',message,'TooltipString',message,'ForegroundColor',color);

function UpdateWarning(h,s,image)
global Stack
u=getappdata(h,'UI'); set(u.Warning,'String','','TooltipString','');
[ready,~]=BackgroundReady(u,s,u.Channel);
if ~ready || ~s.Enabled, return; end
if nargin<3
    if u.Channel>numel(Stack) || isempty(Stack{u.Channel}), return; end
    image=Stack{u.Channel}(:,:,1);
end
try
    report=fSpatialBackground('check',image,s);
    if report.Warn
        message=sprintf('Background mismatch: correlation %.2f; residual variation %.0f%%. Consider another map or lower strength.', ...
            report.Correlation,100*report.ResidualFraction);
        if ~isfinite(report.Correlation), message=report.Message; end
        set(u.Warning,'String',message,'TooltipString',report.Message);
    end
catch err
    set(u.Warning,'String',err.message,'TooltipString',err.message);
end

function Image=ThresholdPreview(Image,Config,hMainGui,channels,frames,Drift)
%FSPATIALBACKGROUNDPREVIEW Threshold native pixels before display transforms.
% Analysis also thresholds each native channel before applying quickwarp.
[filter,background]=strtok(Config.Threshold.Filter,'+');
apply=(strcmp(get(hMainGui.Menu.mApplyCorrections,'Checked'),'on') && ...
    strcmp(get(hMainGui.Menu.mApplyCorrections,'Enable'),'on')) || ...
    (~isempty(hMainGui.Values.PostSpecial) && strcmp(get(hMainGui.ToolBar.ToolChannels(5),'State'),'off'));
for n=1:numel(channels)
    channel=channels(n);
    p=struct('scale',Config.PixSize,'fwhm_estimate',Config.Threshold.FWHM/Config.PixSize, ...
        'binary_image_processing',filter,'background_filter',background);
    p.spatial_background=fSpatialBackground('get',Config.Threshold,channel);
    if strcmp(Config.Threshold.Mode,'relative')
        p.threshold=hMainGui.Values.RelThresh(channel)*1i;
    elseif strcmp(Config.Threshold.Mode,'constant')
        p.threshold=hMainGui.Values.Thresh(channel);
    end
    mask=double(Image2Binary(Image(:,:,n),p));
    frame=frames(min(n,numel(frames)));
    if apply && frame>0 && isreal(frame) && ~isempty(Drift{channel})
        idx=min(frame,size(Drift{channel},3));
        mask=quickwarp(mask,Drift{channel}(:,:,idx),0);
        active=hMainGui.Values.FrameIdx(1);
        if active>1 && ~strcmp(get(hMainGui.Menu.mApplyCorrections,'Checked'),'on') && ...
                ~isempty(hMainGui.Values.PostSpecial)
            idx=min(frame,size(Drift{active},3));
            mask=quickwarp(mask,Drift{active}(:,:,idx),1);
        end
    end
    Image(:,:,n)=mask;
end

function reports=CheckWarnings(threshold,Stack,showDialog)
%FSPATIALBACKGROUNDWARNINGS Check each job once, never on every frame.
if nargin<3, showDialog=true; end
reports=cell(1,numel(Stack));
for channel=1:numel(Stack)
    s=fSpatialBackground('get',threshold,channel);
    if s.Enabled
        reports{channel}=fSpatialBackground('check',Stack{channel}(:,:,1),s);
        if reports{channel}.Warn
            message=sprintf('Channel %d: %s',channel,reports{channel}.Message);
            if showDialog
                fMsgDlg(message,'warning');
            else
                warning('FIESTA:SpatialBackgroundMismatch','%s',message);
            end
        end
    end
end
