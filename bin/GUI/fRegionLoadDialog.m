function [choice,alignRegions] = fRegionLoadDialog
%FREGIONLOADDIALOG Region source selection with opt-in movie alignment.
choice='Cancel';
alignRegions=false;
h=dialog('Name','Load Regions','Tag','FIESTARegionLoadDialog', ...
    'Units','pixels','Position',[100 100 590 175],'Resize','off', ...
    'CloseRequestFcn',@(src,event) uiresume(src));
cleanup=onCleanup(@() DeleteDialog(h));
movegui(h,'center');
setappdata(h,'Choice','Cancel');
uicontrol(h,'Style','text','Units','normalized','Position',[.04 .64 .92 .22], ...
    'String','Load saved regions or convert a folder of ImageJ line ROIs?', ...
    'HorizontalAlignment','left','FontSize',11);
labels={'MAT file','ImageJ ROI folder','Cancel'};
for k=1:numel(labels)
    uicontrol(h,'Style','pushbutton','Units','normalized', ...
        'Position',[.04+(k-1)*.32 .35 .28 .24],'String',labels{k}, ...
        'FontSize',11,'Callback',@(src,event) SelectSource(src));
end
checkbox=uicontrol(h,'Style','checkbox','Units','normalized', ...
    'Position',[.04 .08 .92 .18],'Value',0,'FontSize',11, ...
    'Tag','AlignRegionsToProjection', ...
    'String','align regions with z-projection of current image', ...
    'TooltipString','Apply one shared x/y translation using the maximum z-projection of the active movie channel');
uiwait(h);
if isgraphics(h)
    choice=getappdata(h,'Choice');
    alignRegions=logical(get(checkbox,'Value'));
end

function SelectSource(button)
h=ancestor(button,'figure');
setappdata(h,'Choice',get(button,'String'));
uiresume(h);

function DeleteDialog(h)
if isgraphics(h), delete(h); end
