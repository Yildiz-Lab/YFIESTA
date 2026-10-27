function options=fMotilityStatsDlg(nRegions)
% Compact options shared by selected-scan and pooled-region analysis.
if nargin<1,nRegions=0;end
options=[];c=[50 50 50]/255;fg=[.9 .9 .9];editColor=[.12 .12 .12];
h=dialog('Name','Yildiz Motility Parameters','Tag','hMotilityStatsDlg','Color',c, ...
    'Units','pixels','Position',[100 100 440 520],'Visible','off');movegui(h,'center');
defaults=fYildizMotility('defaults');
label('Thresholds',.93,true);
d=number('Moving distance (pixels)',defaults.dthresh,.865,'dthresh');
t=number('Minimum track span (frames)',defaults.fthresh,.80,'fthresh');
label('Calculate using moving molecules only',.725,true);
v=check('Velocity',defaults.velocity_moverstoggle,.665,'velocity_moverstoggle');
r=check('Run length / interaction time',defaults.runlength_moverstoggle,.61,'runlength_moverstoggle');
l=check('Landing rate',defaults.landingrate_moverstoggle,.555,'landingrate_moverstoggle');
edge=check('Exclude tracks intersecting the dominant end',false,.475,'runlength_countedges');
set(edge,'TooltipString','Exclude velocity, run length and interaction time if any localization lies in the end zone. Landing counts are retained. A direction tie leaves tracks included.');
e=number('End-zone width (pixels)',defaults.edgeregion,.42,'edgeregion');
set(e,'Enable','off');set(edge,'Callback',@(s,ev)set(e,'Enable',onoff(get(s,'Value'))));
pool=check('Pool regions and export automatically',false,.345,'pool');
ids=number(sprintf('Regions (all or 1:5; %d loaded)',nRegions),'all',.285,'regions');
kymo=check('Export kymographs',false,.215,'exportKymographs');
label('Dominant kymograph direction',.155,false);
direction=uicontrol(h,'Style','popupmenu','Units','normalized','Position',[.60 .155 .36 .048], ...
    'String',{'+1 (right)','-1 (left)'},'Tag','kymoDirection','BackgroundColor',editColor,'ForegroundColor',fg,'FontSize',10);
overlay=check('Include molecule overlay',false,.095,'kymoOverlay');
set(direction,'TooltipString','Mirror the image and tracks together to give the chosen dominant direction. Unresolved regions remain unflipped.');
set(overlay,'TooltipString','Save one color TIFF with trajectories drawn in cyan and end-excluded tracks in orange. Uncheck to save the raw grayscale TIFF instead.');
set(pool,'Callback',@togglePool);set(kymo,'Callback',@togglePool);togglePool;
if nRegions==0,set(pool,'Enable','off');end
uicontrol(h,'Style','pushbutton','Units','normalized','Position',[.52 .025 .22 .065], ...
    'String','Run','Tag','Run','BackgroundColor',[.16 .16 .16],'ForegroundColor',fg,'Callback',@accept);
uicontrol(h,'Style','pushbutton','Units','normalized','Position',[.76 .025 .20 .065], ...
    'String','Cancel','BackgroundColor',[.16 .16 .16],'ForegroundColor',fg,'Callback',@(s,ev)delete(h));
set(h,'Visible','on');uiwait(h);if isgraphics(h),delete(h);end

    function label(text,y,bold)
        weight='normal';if bold,weight='bold';end
        uicontrol(h,'Style','text','Units','normalized','Position',[.04 y .69 .05], ...
            'String',text,'HorizontalAlignment','left','FontSize',10,'FontWeight',weight,'BackgroundColor',c,'ForegroundColor',fg);
    end
    function control=number(text,value,y,tag)
        label(text,y,false);if isnumeric(value),value=num2str(value);end
        control=uicontrol(h,'Style','edit','Units','normalized','Position',[.74 y .22 .055], ...
            'String',value,'Tag',tag,'BackgroundColor',editColor,'ForegroundColor',fg,'FontSize',10);
    end
    function control=check(text,value,y,tag)
        control=uicontrol(h,'Style','checkbox','Units','normalized','Position',[.04 y .92 .055], ...
            'String',text,'Value',value,'Tag',tag,'FontSize',10,'BackgroundColor',c,'ForegroundColor',fg);
    end
    function togglePool(varargin)
        enabled=logical(get(pool,'Value'));set(ids,'Enable',onoff(enabled));set(kymo,'Enable',onoff(enabled));
        set([direction overlay],'Enable',onoff(enabled&&get(kymo,'Value')));
    end
    function accept(~,~)
        result=defaults;result.dthresh=str2double(get(d,'String'));result.fthresh=str2double(get(t,'String'));
        result.edgeregion=str2double(get(e,'String'));
        if any(~isfinite([result.dthresh result.fthresh result.edgeregion]))||any([result.dthresh result.fthresh result.edgeregion]<0)
            errordlg('Thresholds must be finite nonnegative numbers.','Motility options');return;
        end
        result.velocity_moverstoggle=logical(get(v,'Value'));result.runlength_moverstoggle=logical(get(r,'Value'));
        result.landingrate_moverstoggle=logical(get(l,'Value'));result.runlength_countedges=logical(get(edge,'Value'));
        result.pool=logical(get(pool,'Value'));result.regions=strtrim(get(ids,'String'));
        result.exportKymographs=result.pool&&logical(get(kymo,'Value'));
        result.kymoDirection=3-2*get(direction,'Value');result.kymoOverlay=logical(get(overlay,'Value'));
        options=result;uiresume(h);
    end
end

function value=onoff(flag)
value='off';if flag,value='on';end
end
