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
global Stack;
global Filament;

fprintf("We are running MTIMBS \n")

FilSelect = [Filament.Selected];
if all(FilSelect==0)
    fMsgDlg('No filaments selected!','error');
    return;
end
MTIMBSFil = Filament(FilSelect==1);

% make an array to store data
% thinking number, length, intensity-background mean value,
% intensity-background std, integral intensity, mean intensity, standard deviation,
% mean background, background standard deviation, ?use this for analysis?
MTIMBSData = nan(length(MTIMBSFil),10);
% which channel would you like to use? Set it to the current one user is on
idx = [real(getFrameIdx(hMainGui)),NaN];

% generate convex boundary list (regions) to check whether points are in it
% or not
if isfield(hMainGui,'Regions') % if no regions, then the whole frame
    m = size(Stack{idx(1)},1); n = size(Stack{idx(1)},2);
    Xbd{1} = [1,m,m,1,1]; Ybd{1} = [1,1,n,n,1];
elseif length(hMainGui.Region) < 1 % if no regions, then the whole frame
    m = size(Stack{idx(1)},1); n = size(Stack{idx(1)},2);
    Xbd{1} = [1,m,m,1,1]; Ybd{1} = [1,1,n,n,1];            
else
    Xbd = cell(1,length(hMainGui.Region)); Ybd = cell(1,length(hMainGui.Region));
    for r=1:length(hMainGui.Region)
        Xbd{r} = hMainGui.Region(r).X; Ybd{r} = hMainGui.Region(r).Y;
    end
end

% the Data store has:
%   x,y coordinates are essential
%   even has intensity measurements at those spots
% MTIMBSFil Data [x, y, z, distance, ?, intensity, ?]
for i=1:length(MTIMBSFil)
    for j = 1:size(MTIMBSFil(i).Results,1) % takes care of different frames for compiled movies
        FilData = MTIMBSFil(i).Data{j};
        X = FilData(:,1)/hMainGui.Values.PixSize; %want in stack coordinates
        Y = FilData(:,2)/hMainGui.Values.PixSize; %want in stack coordinates
        
        % is the shifted line in any region? Will be problematic if regions highly
        % overlap, but should estimate relatively well
        fil_in_region = 0;
        for r=1:length(Xbd)
            pts_in_region = inpolygon(X,Y,Xbd{r},Ybd{r});
            if sum(pts_in_region)/length(pts_in_region) > 0.9 % a certain percentage should be in
                fil_in_region = 1;
                X = X(pts_in_region); Y = Y(pts_in_region); % make sure to only use the pts in the region
            end
        end
        
        if fil_in_region
        MTIMBSData(i,1) = str2double(MTIMBSFil(i).Name(10:end));
        % interpolate like it is a measure (taken from fMainGui "%check for
        % double click"
        LenArea = 0;
        XI=[];
        YI=[];
        for ii=1:length(X)-1
            len=norm([X(ii+1)-X(ii) Y(ii+1)-Y(ii)]);
            LenArea=LenArea+hMainGui.Values.PixSize/1000*norm([X(ii+1)-X(ii) Y(ii+1)-Y(ii)]);
            XI=[XI linspace(X(ii),X(ii+1),ceil(len))];
            YI=[YI linspace(Y(ii),Y(ii+1),ceil(len))];
        end
        %this is how an interpolation for measure is produced
        idx(2) = MTIMBSFil(i).Results(j,1);
        ZI = interp2(double(Stack{idx(1)}(:,:,idx(2))),XI,YI);

        MTIMBSData(i,2) = LenArea;
        MTIMBSData(i,5) = sum(ZI); MTIMBSData(i,6) = mean(ZI); MTIMBSData(i,7) = std(ZI);
        
        %% Do Background Shift and find best one
        % Now we need to do our check background spots and taking the minimum
        % background so just shift X,Y by some set number of pixels and then do
        % it all again
        sp = 10; % how many pixels in any direction to move
        % this table directs how to move the line for background measurement
        % Change mod if you want more or less background measurement comparisons
        mod = [sp 0; -sp 0; 0 sp; 0 -sp; sp sp; sp -sp; -sp sp; -sp -sp];
        
        MTIMBSData(i,8) = 1e6; %high to be quickly updated
        for k = 1:size(mod,1)
            % make the shifted line
            xmod = X+mod(k,1);
            ymod = Y+mod(k,2);
            
            % is the shifted line in any region? Will be problematic if regions highly
            % overlap, but should estimate relatively well
            inregion = 0;
            for r=1:length(Xbd)
                pts_in_region = inpolygon(xmod,ymod,Xbd{r},Ybd{r});
                if sum(pts_in_region)/length(pts_in_region) > 0.9 % a certain percentage should be in
                    inregion = 1;
                    xmod = xmod(pts_in_region); ymod = ymod(pts_in_region); % make sure to only use the pts in the region
                end
            end

            % get intensities only if majorly within the region
            if inregion
                % interpolate like it is a measure (taken from fMainGui "%check for
                % double click"
                LenArea = 0;
                XImod=[];
                YImod=[];
                for ii=1:length(xmod)-1
                    len=norm([X(ii+1)-X(ii) Y(ii+1)-Y(ii)]);
                    XImod=[XImod linspace(xmod(ii),xmod(ii+1),ceil(len))];
                    YImod=[YImod linspace(ymod(ii),ymod(ii+1),ceil(len))];
                end
                %this is how an interpolation for measure is produced
                ZI = interp2(double(Stack{idx(1)}(:,:,idx(2))),XImod,YImod);
                if mean(ZI) < MTIMBSData(i,8)
                    MTIMBSData(i,8) = mean(ZI); MTIMBSData(i,9) = std(ZI);
                end
            end
        end
        end
    end

    
end

% then do intensity - background subtractions (idk how to do standard
% deviation subtraction but that is probably sketch)
MTIMBSData(:,3) = MTIMBSData(:,6)-MTIMBSData(:,8)

fprintf(strcat("Number of MTs: ", num2str(sum(~isnan(MTIMBSData(:,1)))), "\n"))
fprintf(strcat("Microtubule Mean: ", num2str(round(mean(MTIMBSData(:,3),'omitnan'),1)), "\n"))
fprintf(strcat("Length Weight Averaged Mean: ", num2str(round(sum(MTIMBSData(:,3).*MTIMBSData(:,2),'omitnan')/sum(MTIMBSData(:,2),'omitnan'),1)), "\n"))
fprintf(strcat("Microtubule Std: ", num2str(round(std(MTIMBSData(:,3),'omitnan'),1)), "\n"))

% should export a CSV or show a GUI or something