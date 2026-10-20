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
main=getappdata(0,'hMainGui');
count=max(1,numel(Stack));
channel=min(count,main.Values.FrameIdx(1));
h=dialog('Name','Spatial background for thresholding','Tag','SpatialBackgroundDialog', ...
    'Units','pixels','Position',[100 100 640 315],'Resize','off','WindowStyle','normal');
movegui(h,'center');
c=get(h,'Color');
ui.Parent=parent; ui.Threshold=getappdata(parent.fig,'SpatialBackgroundSettings'); ui.Channel=channel;
labels=arrayfun(@(n) sprintf('Channel %d',n),1:count,'UniformOutput',false);
ui.ChannelBox=uicontrol(h,'Style','popupmenu','Units','normalized','Position',[.04 .85 .24 .09], ...
    'String',labels,'Value',channel,'Tag','SpatialBackgroundChannel','Callback',@(src,event) ChangeChannel(h));
ui.Enabled=uicontrol(h,'Style','checkbox','Units','normalized','Position',[.32 .85 .64 .09], ...
    'String','Use spatial background image','Tag','SpatialBackgroundEnabled','BackgroundColor',c);
ui.File=uicontrol(h,'Style','text','Units','normalized','Position',[.04 .67 .92 .13], ...
    'HorizontalAlignment','left','BackgroundColor',c,'Tag','SpatialBackgroundFile');
uicontrol(h,'Style','pushbutton','Units','normalized','Position',[.04 .53 .26 .10], ...
    'String','Load background TIFF...','Tag','SpatialBackgroundLoad','Callback',@(src,event) LoadMap(h));
uicontrol(h,'Style','pushbutton','Units','normalized','Position',[.32 .53 .16 .10], ...
    'String','Clear','Callback',@(src,event) ClearMap(h));
uicontrol(h,'Style','text','Units','normalized','Position',[.04 .36 .23 .10], ...
    'String','Camera baseline (counts)','BackgroundColor',c,'HorizontalAlignment','left');
ui.Baseline=uicontrol(h,'Style','edit','Units','normalized','Position',[.28 .37 .14 .09], ...
    'Tag','SpatialBackgroundBaseline');
uicontrol(h,'Style','text','Units','normalized','Position',[.49 .36 .29 .10], ...
    'String','Correction strength (1 = full)','BackgroundColor',c,'HorizontalAlignment','left');
ui.Strength=uicontrol(h,'Style','edit','Units','normalized','Position',[.80 .37 .14 .09], ...
    'Tag','SpatialBackgroundStrength');
uicontrol(h,'Style','text','Units','normalized','Position',[.04 .17 .92 .16], ...
    'String',{'1 = median reference signal above the baseline.','The pattern scales to each movie frame; raw pixels remain unchanged for fitting.'}, ...
    'BackgroundColor',c,'HorizontalAlignment','left');
uicontrol(h,'Style','pushbutton','Units','normalized','Position',[.04 .04 .26 .10], ...
    'String','Preview current frame','Tag','SpatialBackgroundPreview','Callback',@(src,event) Preview(h));
uicontrol(h,'Style','pushbutton','Units','normalized','Position',[.58 .04 .17 .10], ...
    'String','Apply','Tag','SpatialBackgroundApply','Callback',@(src,event) Apply(h));
uicontrol(h,'Style','pushbutton','Units','normalized','Position',[.78 .04 .17 .10], ...
    'String','Cancel','Callback',@(src,event) delete(h));
setappdata(h,'UI',ui); Populate(h);

function Populate(h)
u=getappdata(h,'UI'); s=fSpatialBackground('get',u.Threshold,u.Channel);
set(u.Enabled,'Value',s.Enabled); set(u.Baseline,'String',num2str(s.Baseline)); set(u.Strength,'String',num2str(s.Strength));
if isempty(s.NormalizedMap)
    text='No background image loaded.';
else
    [~,name,ext]=fileparts(s.File);
    text=sprintf('%s%s | %d x %d | %d TIFF page(s) averaged | reference median %.2f counts', ...
        name,ext,size(s.NormalizedMap,2),size(s.NormalizedMap,1),s.FrameCount,s.SourceMedian);
end
set(u.File,'String',text,'TooltipString',s.File);

function s=Collect(h)
u=getappdata(h,'UI'); s=fSpatialBackground('get',u.Threshold,u.Channel);
baseline=str2double(get(u.Baseline,'String')); strength=str2double(get(u.Strength,'String'));
validateattributes(baseline,{'numeric'},{'real','scalar','finite','nonnegative'});
validateattributes(strength,{'numeric'},{'real','scalar','finite','nonnegative'});
if baseline~=s.Baseline && ~isempty(s.NormalizedMap)
    s=fSpatialBackground('load',s.File,baseline,strength);
end
s.Baseline=baseline; s.Strength=strength; s.Enabled=logical(get(u.Enabled,'Value'));

function ChangeChannel(h)
u=getappdata(h,'UI');
try
    s=Collect(h); u.Threshold=fSpatialBackground('set',u.Threshold,u.Channel,s);
    u.Channel=get(u.ChannelBox,'Value'); setappdata(h,'UI',u); Populate(h);
catch err
    set(u.ChannelBox,'Value',u.Channel); fMsgDlg(err.message,'error');
end

function LoadMap(h)
global Stack
u=getappdata(h,'UI');
[file,folder]=uigetfile({'*.tif;*.tiff','Background TIFF (*.tif, *.tiff)'},'Load Gaussian-blurred background image');
if isequal(file,0), return; end
try
    s=Collect(h);
    s=fSpatialBackground('load',fullfile(folder,file),s.Baseline,s.Strength);
    if u.Channel<=numel(Stack)
        fSpatialBackground('validate',s,[size(Stack{u.Channel},1) size(Stack{u.Channel},2)],'constant');
    end
    u.Threshold=fSpatialBackground('set',u.Threshold,u.Channel,s);
    setappdata(h,'UI',u); Populate(h);
catch err
    fMsgDlg(err.message,'error');
end

function ClearMap(h)
u=getappdata(h,'UI');
u.Threshold=fSpatialBackground('set',u.Threshold,u.Channel,fSpatialBackground('defaults'));
setappdata(h,'UI',u); Populate(h);

function Apply(h)
global Stack
u=getappdata(h,'UI');
try
    s=Collect(h);
    mode='constant';
    if get(u.Parent.pGeneral.rVariable,'Value'), mode='variable'; end
    imageSize=[];
    if u.Channel<=numel(Stack), imageSize=[size(Stack{u.Channel},1) size(Stack{u.Channel},2)]; end
    fSpatialBackground('validate',s,imageSize,mode);
    u.Threshold=fSpatialBackground('set',u.Threshold,u.Channel,s);
    setappdata(u.Parent.fig,'SpatialBackgroundSettings',u.Threshold);
    if fSpatialBackground('enabled',u.Threshold,1:max(1,numel(Stack)))
        set(u.Parent.pGeneral.cFilterBackground,'Value',0);
    end
    if s.Enabled && u.Channel<=numel(Stack)
        report=fSpatialBackground('check',Stack{u.Channel}(:,:,1),s);
        if report.Warn, fMsgDlg(report.Message,'warning'); end
    end
    delete(h);
catch err
    fMsgDlg(err.message,'error');
end

function Preview(h)
global Stack
u=getappdata(h,'UI');
try
    s=Collect(h);
    main=getappdata(0,'hMainGui'); channel=u.Channel;
    if channel>numel(Stack) || isempty(Stack{channel})
        error('FIESTA:NoImage','Load a movie before previewing the background correction.');
    end
    mode='constant'; if get(u.Parent.pGeneral.rRelative,'Value'), mode='relative'; end
    if get(u.Parent.pGeneral.rVariable,'Value'), mode='variable'; end
    test=s; test.Enabled=true;
    fSpatialBackground('validate',test,[size(Stack{channel},1) size(Stack{channel},2)],mode);
    if ~s.Enabled
        error('FIESTA:SpatialBackgroundDisabled','Enable Use spatial background image to preview this correction.');
    end
    frame=main.Values.FrameIdx(min(channel+1,numel(main.Values.FrameIdx)));
    if ~isreal(frame) || frame<1, frame=1; end
    frame=min(frame,size(Stack{channel},3)); image=double(Stack{channel}(:,:,frame));
    filter='none';
    if get(u.Parent.pGeneral.rFilterAverage,'Value'), filter='average'; end
    if get(u.Parent.pGeneral.rFilterSmooth,'Value'), filter='smooth'; end
    p=struct('spatial_background',test,'binary_image_processing',filter,'background_filter','', ...
        'threshold',main.Values.Thresh(channel));
    if strcmp(mode,'relative'), p.threshold=main.Values.RelThresh(channel)*1i; end
    [mask,details]=fSpatialBackground('threshold',image,p);
    preview=figure('Name',sprintf('Spatial background preview - channel %d, frame %d',channel,frame), ...
        'NumberTitle','off','Tag','SpatialBackgroundPreviewFigure');
    subplot(2,2,1); imagesc(image); axis image; colorbar; title('Raw frame (counts)');
    subplot(2,2,2); imagesc(s.NormalizedMap); axis image; colorbar; title('Normalized background (median = 1)');
    subplot(2,2,3); imagesc(details.ThresholdMap); axis image; colorbar; title('Pixel threshold (counts)');
    subplot(2,2,4); imagesc(mask,[0 1]); axis image; title(sprintf('Detection mask (strength %.2g)',test.Strength));
    colormap(preview,'gray');
    report=fSpatialBackground('check',image,test);
    if report.Warn
        annotation(preview,'textbox',[.02 .01 .96 .07],'String',['Background mismatch: ' report.Message], ...
            'EdgeColor','none','Color',[.9 .35 .05],'Interpreter','none','FontSize',9);
    end
catch err
    fMsgDlg(err.message,'error');
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
