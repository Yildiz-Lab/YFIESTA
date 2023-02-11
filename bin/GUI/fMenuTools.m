function fMenuTools(func,varargin)
switch func
    case 'MeasureLine'
        MeasureLine(varargin{1});
    case 'MeasureSegLine'
        MeasureSegLine(varargin{1});
    case 'MeasureFreehand'
        MeasureFreehand(varargin{1});
    case 'MeasureRect'
        MeasureRect(varargin{1});       
    case 'MeasureEllipse'
        MeasureEllipse(varargin{1});
    case 'MeasurePolygon'
        MeasurePolygon(varargin{1});
    case 'ScanLine'
        ScanLine(varargin{1});
    case 'ScanSegLine'
        ScanSegLine(varargin{1});
    case 'ScanFreehand'
        ScanFreehand(varargin{1});      
    case 'ShowKymoGraph'
        ShowKymoGraph(varargin{1});   
    case 'MTIMBS'
        MTIMBS(varargin{1});
end

function MeasureLine(hMainGui)
fRightPanel('ToolsMeasurePanel',hMainGui);
fToolBar('Cursor',hMainGui);
hMainGui=getappdata(0,'hMainGui');
set(hMainGui.RightPanel.pTools.bLine,'Value',1);
hMainGui.CursorMode=get(hMainGui.RightPanel.pTools.bLine,'UserData');
hMainGui.Values.CursorDownPos(:)=0;
setappdata(0,'hMainGui',hMainGui);

function MeasureSegLine(hMainGui)
fRightPanel('ToolsMeasurePanel',hMainGui);
fToolBar('Cursor',hMainGui);
hMainGui=getappdata(0,'hMainGui');
set(hMainGui.RightPanel.pTools.bSegLine,'Value',1);
hMainGui.CursorMode=get(hMainGui.RightPanel.pTools.bSegLine,'UserData');
hMainGui.Values.CursorDownPos(:)=0;
setappdata(0,'hMainGui',hMainGui);

function MeasureFreehand(hMainGui)
fRightPanel('ToolsMeasurePanel',hMainGui);
fToolBar('Cursor',hMainGui);
hMainGui=getappdata(0,'hMainGui');
set(hMainGui.RightPanel.pTools.bFreehand,'Value',1);
hMainGui.CursorMode=get(hMainGui.RightPanel.pTools.bFreehand,'UserData');
hMainGui.Values.CursorDownPos(:)=0;
setappdata(0,'hMainGui',hMainGui);

function MeasureRect(hMainGui)
fRightPanel('ToolsMeasurePanel',hMainGui);
fToolBar('Cursor',hMainGui);
hMainGui=getappdata(0,'hMainGui');
set(hMainGui.RightPanel.pTools.bRectangle,'Value',1);
hMainGui.CursorMode=get(hMainGui.RightPanel.pTools.bRectangle,'UserData');
hMainGui.Values.CursorDownPos(:)=0;
setappdata(0,'hMainGui',hMainGui);

function MeasureEllipse(hMainGui)
fRightPanel('ToolsMeasurePanel',hMainGui);
fToolBar('Cursor',hMainGui);
hMainGui=getappdata(0,'hMainGui');
set(hMainGui.RightPanel.pTools.bEllipse,'Value',1);
hMainGui.CursorMode=get(hMainGui.RightPanel.pTools.bEllipse,'UserData');
hMainGui.Values.CursorDownPos(:)=0;
setappdata(0,'hMainGui',hMainGui);

function MeasurePolygon(hMainGui)
fRightPanel('ToolsMeasurePanel',hMainGui);
fToolBar('Cursor',hMainGui);
hMainGui=getappdata(0,'hMainGui');
set(hMainGui.RightPanel.pTools.bPolygon,'Value',1);
hMainGui.CursorMode=get(hMainGui.RightPanel.pTools.bPolygon,'UserData');
hMainGui.Values.CursorDownPos(:)=0;
setappdata(0,'hMainGui',hMainGui);

function ScanLine(hMainGui)
fRightPanel('ToolsScanPanel',hMainGui);
fToolBar('Cursor',hMainGui);
hMainGui=getappdata(0,'hMainGui');
set(hMainGui.RightPanel.pTools.bLineScan,'Value',1);
hMainGui.CursorMode=get(hMainGui.RightPanel.pTools.bLineScan,'UserData');
hMainGui.Values.CursorDownPos(:)=0;
setappdata(0,'hMainGui',hMainGui);

function ScanSegLine(hMainGui)
fRightPanel('ToolsScanPanel',hMainGui);
fToolBar('Cursor',hMainGui);
hMainGui=getappdata(0,'hMainGui');
set(hMainGui.RightPanel.pTools.bSegLineScan,'Value',1);
hMainGui.CursorMode=get(hMainGui.RightPanel.pTools.bSegLineScan,'UserData');
hMainGui.Values.CursorDownPos(:)=0;
setappdata(0,'hMainGui',hMainGui);

function ScanFreehand(hMainGui)
fRightPanel('ToolsScanPanel',hMainGui);
fToolBar('Cursor',hMainGui);
hMainGui=getappdata(0,'hMainGui');
hMainGui.CursorMode='normal';
fRightPanel('CreateFilamentScan',hMainGui);
hMainGui=getappdata(0,'hMainGui');
hMainGui.Values.CursorDownPos(:)=0;
setappdata(0,'hMainGui',hMainGui);

function ShowKymoGraph(hMainGui)
fRightPanel('ToolsScanPanel',hMainGui);
fToolBar('Cursor',hMainGui);
hMainGui=getappdata(0,'hMainGui');
set(hMainGui.RightPanel.pTools.cShowKymoGraph,'Value',1);
fRightPanel('ShowKymoGraph',hMainGui);
hMainGui.Values.CursorDownPos(:)=0;
setappdata(0,'hMainGui',hMainGui);

% JS Edit 2023/02/10 get MTIMBS in FIESTA
function MTIMBS(hMainGui)
global Config;
global Stack;
global Filament;

fprintf("We are running MTIMBS \n")

FilSelect = [Filament.Selected];
if all(FilSelect==0)
    fMsgDlg('No filaments selected!','error');
    return;
end
MTIMBSFil = Filament(FilSelect==1);

% the Data store has:
%   x,y coordinates are essential
%   even has intensity measurements at those spots
% MTIMBSFil Data [x, y, z, distance, ?, intensity, ?]
for i=1:length(MTIMBSFil)
    size(MTIMBSFil(i).Data,2)
    FilData = MTIMBSFil(i).Data{1};
    X = FilData(:,1);
    Y = FilData(:,2);
    % interpolate for every pixel
    t = 1:length(X);
    lint = linspace(1,length(X),round(4*max(FilData(:,4))/Config.PixSize,0));
    xchk = interp1(t, X, lint)/Config.PixSize;     ychk = interp1(t, Y, lint)/Config.PixSize;
    TXY = unique(horzcat(transpose(round(xchk,0)), transpose(round(ychk,0))),'rows','stable');
    
    sumI = 0;
    for j=1:size(TXY,1)
    sumI = sumI + double(Stack{2}(TXY(j,1),TXY(j,2)))
    end
    sumI/size(TXY,1)
    mean(FilData(:,6),'omitnan')
end


