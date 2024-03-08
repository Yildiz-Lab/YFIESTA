function fPathStatsGui(func,varargin)
switch func
    case 'Create'
        Create;
    case 'Update'
        Update(varargin{1});        
    case 'Draw'
        Draw(varargin{1});  
    case 'bToggleToolCursor'
        bToggleToolCursor(varargin{1});  
    case 'bToolPan'
        bToolPan(varargin{1});
    case 'bToolZoomIn'
        bToolZoomIn(varargin{1});
    case 'Drift'
        Drift(varargin{1});  
    case 'Load'
        Load(varargin{1});  
    case 'Cancel'
        Cancel(varargin{1});          
    case 'Ok'
        Ok(varargin{1}); 
    case 'Crop'
        Crop(varargin{1});
    case 'FIONA'
        FIONA(varargin{1});
    case 'StepStats'
        StepStats(varargin{1});
end

function Create
global Molecule;
global Filament;
global Index;
global Config;
global FiestaDir;

% JS Edit for dark theme
ctext = 0.9*ones(1,3);
cbutton = [0.16,0.16,0.16];

h=findobj('Tag','hPathsStatsGui');
close(h)

MolSelect = [Molecule.Selected];
FilSelect = [Filament.Selected];
if all(MolSelect==0) && all(FilSelect==0)
    fMsgDlg('No track selected!','error');
    return;
end
options = fPathStatsDlg;
if isempty(options)
    return;
end
if strcmp(options.mode,'average')
    AverageDis = options.dis;
else
    AverageDis  = 0;
end

if strcmp(options.mode,'filament')
    if isempty(Filament)
        fMsgDlg('No filaments present!','error');
        return;
    end
    if options.align == 1
        fMenuStatistics('AlignFilament');
    end
    if all(FilSelect==0)
        PathFilSelect = ones(size(Filament));
    else
        PathFilSelect = zeros(size(Filament));
        for n=1:length(Filament)
            if size(Filament(n).Results,1)==1
                PathFilSelect(n)=1;
            end
        end
        if all(PathFilSelect==1)
            FilSelect(:)=0;
        end
    end
    if all(MolSelect==0) && all(PathFilSelect==1)
        fMsgDlg('No molecules selected!','error');
        return;
    end
    PathFilSelect = find(PathFilSelect==1);
else
    PathFilSelect =[];
end
PathMol = Molecule(MolSelect==1);
PathMol = rmfield(PathMol,'Type');
PathFil = Filament(FilSelect==1);
PathFil = rmfield(PathFil,{'PosStart','PosCenter','PosEnd','Data'});
PathStats = [PathMol PathFil];

Index = [ find(MolSelect==1) find(FilSelect==1)*1i ];

hPathsStatsGui.fig = figure('Units','normalized','DockControls','off','IntegerHandle','off','MenuBar','none','Name','Path Statistics',...
                      'NumberTitle','off','HandleVisibility','callback','Tag','hPathsStatsGui',...
                      'Visible','off','Resize','off','WindowStyle','modal');
                  
if ispc || ismac
    set(hPathsStatsGui.fig,'Color',[50 50 50]/255);
end

c=get(hPathsStatsGui.fig,'Color');

hPathsStatsGui.pPlotXYZPanel = uipanel('Parent',hPathsStatsGui.fig,'Position',[0.45 0.625 0.525 0.35],'Tag','PlotPanel','BackgroundColor',c);

hPathsStatsGui.tCalcPath = uicontrol('Parent',hPathsStatsGui.pPlotXYZPanel,'Units','normalized','Position',[00 0.3 1 0.4],'FontSize',12,'FontWeight','bold',...
                        'String','Calculating Path','Style','text','Tag','tCalcPath','HorizontalAlignment','center','BackgroundColor',c,'ForegroundColor',ctext);
                    
hPathsStatsGui.aPlotXYZ = axes('Parent',hPathsStatsGui.pPlotXYZPanel,'Units','normalized','OuterPosition',[0 0 1 1],'Tag','aPlotXYZ','Visible','off');

% JS Edit 2022/06/03 for step sorting
hPathsStatsGui.pMAll = uipanel('Parent',hPathsStatsGui.fig,'Units','normalized','Title','List',...
                             'Position',[0.025 0.83 0.375 0.15],'Tag','pMAll','BackgroundColor',c,'ForegroundColor',ctext);

% hPathsStatsGui.lAll = uicontrol('Parent',hPathsStatsGui.fig,'Units','normalized','BackgroundColor',[1 1 1],'Callback','fPathStatsGui(''Update'',getappdata(0,''hPathsStatsGui''));',...
%                           'Position',[0.025 0.83 0.375 0.15],'String',{PathStats.Name},'Style','listbox','Value',1,'Tag','lAll','Max',length(PathStats));                         

hPathsStatsGui.lAll = uicontrol('Parent',hPathsStatsGui.pMAll,'Units','normalized','BackgroundColor',[1 1 1],'Callback','fPathStatsGui(''Update'',getappdata(0,''hPathsStatsGui''));',...
                           'Position',[0.025 0.05 0.8 0.9],'String',{PathStats.Name},'Style','listbox','Value',1,'Tag','lAll','Max',length(PathStats),'BackgroundColor',c,'ForegroundColor',ctext);                         

hPathsStatsGui.bStepSort = uicontrol('Parent',hPathsStatsGui.pMAll,'Units','normalized','Callback','fPathStatsGui(''Update'',getappdata(0,''hPathsStatsGui''));',...
                            'Position',[0.8375 0.3 0.15 0.4],'String','Sort Best','Tag','bReset','BackgroundColor',cbutton,'ForegroundColor',ctext); 
% end of JS Edit 2022/06/03 for step sorting
                        
hPathsStatsGui.pOptions = uipanel('Parent',hPathsStatsGui.fig,'Units','normalized','Title','Options',...
                             'Position',[0.025 0.625 0.375 0.19],'Tag','pOptions','BackgroundColor',c,'ForegroundColor',ctext);
                        
hPathsStatsGui.bAuto = uicontrol('Parent',hPathsStatsGui.pOptions,'Units','normalized','Callback','fPathStatsGui(''Update'',getappdata(0,''hPathsStatsGui''));',...
                            'Position',[0.75 0.79 0.225 0.23],'String','Auto Fit','Tag','bReset','BackgroundColor',cbutton,'ForegroundColor',ctext);    
                        
hPathsStatsGui.bReset = uicontrol('Parent',hPathsStatsGui.pOptions,'Units','normalized','Callback','fPathStatsGui(''Update'',getappdata(0,''hPathsStatsGui''));',...
                            'Position',[0.75 0.535 0.225 0.23],'String','Reset plots','Tag','bReset','BackgroundColor',cbutton,'ForegroundColor',ctext);                        
                        
hPathsStatsGui.bDisregard = uicontrol('Parent',hPathsStatsGui.pOptions,'Units','normalized','Callback','fPathStatsGui(''Update'',getappdata(0,''hPathsStatsGui''));',...
                            'Position',[0.75 0.28 0.225 0.23],'String','Disregard','Tag','bReset','BackgroundColor',cbutton,'ForegroundColor',ctext);                            
                        
% JS Edit 2022/05/29 to add direct line to FIONA 
hPathsStatsGui.bFiona = uicontrol('Parent',hPathsStatsGui.pOptions,'Units','normalized',...
                            'Callback','fPathStatsGui(''FIONA'',getappdata(0,''hPathsStatsGui''));',...
                            'Position',[0.75 0.025 0.225 0.23],'String','FIONA stepping','Tag','bFIONA','BackgroundColor',cbutton,'ForegroundColor',ctext); 

%JS Edit 2022/09/28 for ability to plot neighbors in other channel
hPathsStatsGui.cNeighbor = uicontrol('Parent',hPathsStatsGui.pOptions,'Units','normalized','Position',[0.52 0.025 0.225 0.23],'Enable','on',...
                                        'String','Show Neighbors','Style','checkbox','Tag','neighbors','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext,'Value',1);  
% End of JS Edit 2022/05/29
                        
hPathsStatsGui.rLinear = uicontrol('Parent',hPathsStatsGui.pOptions,'Units','normalized','Callback','fPathStatsGui(''Update'',getappdata(0,''hPathsStatsGui''));',...
                            'Position',[0.1 0.84 0.6 0.15],'String','Linear path','Style','radiobutton','BackgroundColor',c,'ForegroundColor',ctext,'Tag','rLinear','Value',0);                         

hPathsStatsGui.rPoly2 = uicontrol('Parent',hPathsStatsGui.pOptions,'Units','normalized','Callback','fPathStatsGui(''Update'',getappdata(0,''hPathsStatsGui''));',...
                            'Position',[0.1 0.68 0.6 0.15],'String','2nd deg polynomial path','Style','radiobutton','BackgroundColor',c,'ForegroundColor',ctext,'Tag','rPoly2','Value',0);          

hPathsStatsGui.rPoly3 = uicontrol('Parent',hPathsStatsGui.pOptions,'Units','normalized','Callback','fPathStatsGui(''Update'',getappdata(0,''hPathsStatsGui''));',...
                            'Position',[0.1 0.52 0.6 0.15],'String','3rd deg polynomial path','Style','radiobutton','BackgroundColor',c,'ForegroundColor',ctext,'Tag','rPoly3','Value',0);                          

hPathsStatsGui.rFilament = uicontrol('Parent',hPathsStatsGui.pOptions,'Units','normalized','Callback','fPathStatsGui(''Update'',getappdata(0,''hPathsStatsGui''));',...
                                    'Position',[0.1 0.36 0.6 0.15],'String','Filament path','Style','radiobutton','BackgroundColor',c,'ForegroundColor',ctext,'Tag','rFilament','Value',0);                          

if isempty(Filament)
    set(hPathsStatsGui.rFilament,'Enable','off');
end

hPathsStatsGui.rAverage = uicontrol('Parent',hPathsStatsGui.pOptions,'Units','normalized','Callback','fPathStatsGui(''Update'',getappdata(0,''hPathsStatsGui''));',...
                            'Position',[0.1 0.2 0.6 0.15],'String','Average path','Style','radiobutton','BackgroundColor',c,'ForegroundColor',ctext,'Tag','rAverage','Value',0);   

hPathsStatsGui.tRegion = uicontrol('Parent',hPathsStatsGui.pOptions,'Units','normalized','Enable','off','HorizontalAlignment','left',...
                              'Position',[0.03 0.02 0.15 0.14],'String','Region:','Style','text','Tag','tRegion','BackgroundColor',c,'ForegroundColor',ctext);                         

hPathsStatsGui.eAverage = uicontrol('Parent',hPathsStatsGui.pOptions,'Units','normalized','Callback','fPathStatsGui(''Update'',getappdata(0,''hPathsStatsGui''));','Enable','off',...
                              'Position',[0.12 0.02 0.3 0.14],'String','1000','FontSize',8,'Style','edit','Tag','eAverage','BackgroundColor',[1 1 1]);                         
                          
hPathsStatsGui.tNM = uicontrol('Parent',hPathsStatsGui.pOptions,'Units','normalized','Position',[0.45 0.02 0.05 0.14],'Enable','off',...
                               'String','nm','Style','text','Tag','tNM','HorizontalAlignment','left','BackgroundColor',c);
                           
% JS Edit 2023/01/09 Crop option
hPathsStatsGui.bCrop = uicontrol('Parent',hPathsStatsGui.pOptions,'Units','normalized','Callback','fPathStatsGui(''Crop'',getappdata(0,''hPathsStatsGui''));',...
                            'Position',[0.5 0.28 0.225 0.23],'String','Crop','Tag','bCrop','BackgroundColor',cbutton,'ForegroundColor',ctext);                            
% End of JS Edit 2023/01/09
                        
hPathsStatsGui.pPlotDistPanel = uipanel('Parent',hPathsStatsGui.fig,'Position',[0.025 0.28 0.95 0.335],'Tag','PlotPanel','BackgroundColor',c,'ForegroundColor',ctext);

hPathsStatsGui.aPlotDist = axes('Parent',hPathsStatsGui.pPlotDistPanel,'Units','normalized','OuterPosition',[0 0 1 1],'Tag','aPlotDist');

hPathsStatsGui.pPlotSidePanel = uipanel('Parent',hPathsStatsGui.fig,'Position',[0.025 0.05 0.95 0.22],'Tag','PlotPanel','BackgroundColor',c,'ForegroundColor',ctext);

hPathsStatsGui.aPlotSide = axes('Parent',hPathsStatsGui.pPlotSidePanel ,'Units','normalized','OuterPosition',[0 0 1 1],'Tag','aPlotSide');
    
% JS Edit 2022/06/03 to compile stats
hPathsStatsGui.bPlotStepStats = uicontrol('Parent',hPathsStatsGui.fig,'Units','normalized','Callback','fPathStatsGui(''StepStats'',getappdata(0,''hPathsStatsGui''));',...
                            'Position',[0.1 0.01 0.175 0.03],'String','Compile Stepping Stats','Tag','bStepStats','BackgroundColor',cbutton,'ForegroundColor',ctext);
% end of JS Edit 2022/06/03

hPathsStatsGui.bOk = uicontrol('Parent',hPathsStatsGui.fig,'Units','normalized','Callback','fPathStatsGui(''Ok'',getappdata(0,''hPathsStatsGui''));',...
                            'Position',[0.575 0.01 0.175 0.03],'String','Ok','Tag','bOk','BackgroundColor',cbutton,'ForegroundColor',ctext);

hPathsStatsGui.bCancel = uicontrol('Parent',hPathsStatsGui.fig,'Units','normalized','Callback','fPathStatsGui(''Cancel'',getappdata(0,''hPathsStatsGui''));',...
                             'Position',[0.8 0.01 0.175 0.03],'String','Cancel','Tag','bCancel','BackgroundColor',cbutton,'ForegroundColor',ctext);                        

setappdata(0,'hPathsStatsGui',hPathsStatsGui);
setappdata(hPathsStatsGui.fig,'PathStats',PathStats);
setappdata(hPathsStatsGui.fig,'PathFilSelect',PathFilSelect);

if ~isempty(PathStats)
    if strcmp(options.mode,'filament')
        PosFil = InterpolFil(Filament(PathFilSelect));
    else
        PosFil = [];      
    end    
    setappdata(hPathsStatsGui.fig,'PosFil',PosFil);       
    nPathStats = length(PathStats);
    hMainGui=getappdata(0,'hMainGui');
    dirStatus = [FiestaDir.AppData 'fiestastatus' filesep];  
    parallelprogressdlg('String','Calculating path','Max',nPathStats,'Parent',hMainGui.fig,'Directory',FiestaDir.AppData);
    [PathStats(:).AverageDis] = deal(0);
    fil = Filament(PathFilSelect);
    parfor (n=1:nPathStats,Config.NumCores)
        Path = [];
        if isempty(PathStats(n).PathData) || options.overwrite == 1
            if size(PathStats(n).Results,1)> 9 && ~strcmp(options.mode,'filament')
                if strcmp(options.mode,'fit')
                    if strcmp(options.fit,'auto')
                        [param1,resnorm1] = PathFitLinear(PathStats(n).Results(:,3:5));
                        [param2,resnorm2] = PathFitCurved(PathStats(n).Results(:,3:5),3);
                        [param3,resnorm3] = PathFitCurved(PathStats(n).Results(:,3:5),4);
                        if all(resnorm1*1.1<[resnorm2 resnorm3])
                            Path = EvalLinearPath(param1,PathStats(n).Results(:,3:5));
                            PathStats(n).AverageDis = -1;
                        elseif resnorm2*1.1<resnorm3
                            Path = EvalCurvedPath(param2,PathStats(n).Results(:,3:5));
                            PathStats(n).AverageDis = -2;
                        else
                            Path = EvalCurvedPath(param3,PathStats(n).Results(:,3:5));
                            PathStats(n).AverageDis = -3;
                        end
                    elseif strcmp(options.fit,'poly1')
                        [param1,~] = PathFitLinear(PathStats(n).Results(:,3:5));
                        Path = EvalLinearPath(param1,PathStats(n).Results(:,3:5));
                        PathStats(n).AverageDis = -1;
                    elseif strcmp(options.fit,'poly2')
                        [param2,~] = PathFitCurved(PathStats(n).Results(:,3:5),3); 
                        Path = EvalCurvedPath(param2,PathStats(n).Results(:,3:5));
                        PathStats(n).AverageDis = -2;
                    else
                        [param3,~] = PathFitCurved(PathStats(n).Results(:,3:5),4);
                        Path = EvalCurvedPath(param3,PathStats(n).Results(:,3:5));
                        PathStats(n).AverageDis = -3;
                    end
                elseif strcmp(options.mode,'average')
                    Path = AveragePath(PathStats(n).Results(:,1:5),AverageDis);
                    PathStats(n).AverageDis = AverageDis; 
                end
            elseif strcmp(options.mode,'filament')    
                nFil = size(PosFil,1);
                if nFil > 1
                    s = ones(nFil,1)*Inf;
                    for m = 1:nFil
                        s(m) = GetFilament(double(PathStats(n).Results(:,1:5)),fil(m));
                    end
                    [~,k] = min(s);
                else
                    k = 1;
                end
                Path = EvalFilamentPath(PathStats(n).Results,PosFil(k,:),fil(k));
                if options.refboth == 0 
                    Path = real(Path);
                end
                PathStats(n).AverageDis = -4;
            else
                if size(PathStats(n).Results,1)>2
                    [param1,~] = PathFitLinear(PathStats(n).Results(:,3:5));
                    Path = EvalLinearPath(param1,PathStats(n).Results(:,3:5));
                    PathStats(n).AverageDis = -1;
                else
                    Path = PathStats(n).Results(:,3:6);
                    Path(:,5:6) = NaN;
                    PathStats(n).AverageDis=-1;
                end
            end
            PathStats(n).PathData=Path
        else
            PathStats(n).AverageDis=-5;
        end
        fSave(dirStatus,n);
    end
    parallelprogressdlg('close');
    if PathStats(1).AverageDis == -1
        set(hPathsStatsGui.rLinear,'Value',1);   
    elseif PathStats(1).AverageDis == -2
        set(hPathsStatsGui.rPoly2,'Value',1);   
    elseif PathStats(1).AverageDis == -3
        set(hPathsStatsGui.rPoly3,'Value',1); 
    elseif PathStats(1).AverageDis == -4
        set(hPathsStatsGui.rFilament,'Value',1); 
    elseif PathStats(1).AverageDis > 0
        set(hPathsStatsGui.rAverage,'Value',1); 
        set(hPathsStatsGui.eAverage,'Enable','on','String',num2str(PathStats(1).AverageDis)); 
        set(hPathsStatsGui.tNM,'Enable','on');
    end
    setappdata(hPathsStatsGui.fig,'PathStats',PathStats);
    Draw(hPathsStatsGui);    
end
fPlaceFig(hPathsStatsGui.fig,'big');
  
function path = InterpolFil(Filament)
global Config;
global FiestaDir;
nFil = length(Filament);
hMainGui=getappdata(0,'hMainGui');
nData = zeros(1,nFil);
for n = 1:nFil
    nData(n) = length(Filament(n).Data);
end
dirStatus = [FiestaDir.AppData 'fiestastatus' filesep];  
parallelprogressdlg('String','Interpolating Filaments','Max',nFil,'Parent',hMainGui.fig,'Directory',FiestaDir.AppData);
path = cell(nFil,max(nData));
for n = 1:nFil
    filData = Filament(n).Data;
    parfor (m = 1:nData(n),Config.NumCores)
        X = double(filData{m}(:,1));
        Y = double(filData{m}(:,2));
        Z = double(filData{m}(:,3));
        if numel(X)>1
            P = 1:length(X);
            pi = -5:0.001:length(X)+5;
            pX = interp1(P,X,pi,'linear','extrap');
            pY = interp1(P,Y,pi,'linear','extrap');
            pZ = [];
            if ~any(isnan(Z))
                pZ = interp1(P,Z,pi,'linear','extrap');
            end
        else
            pX = X;
            pY = Y;
            pZ = [];
            if ~isnan(Z)
                pZ = Z;
            end
        end
        path{n,m} = [pX' pY' pZ'];
    end
    fSave(dirStatus,n);
end
parallelprogressdlg('close');

function res = GetFilament(Results,Filament)
if length(Filament.Data) == 1
    if any(isnan(Results(:,5)))
        XY = double(Filament.Data{1}(:,1:2));
        res = min(pdist2(XY,Results(:,3:4)));
    else
        XYZ = double(Filament.Data{1}(:,1:3));
        res = min(pdist2(XYZ,Results(:,3:5)));
    end
else
    nData = size(Results,1);
    res = zeros(1,nData);
    for n = 1:nData
        [~,k] = min(abs(Results(n,1)-Filament.Results(:,1)));
        idx = k(1);
        if any(isnan(Results(:,5)))
            XY = Filament.Data{idx}(:,1:2);
            res(n) = min(pdist2(XY,Results(n,3:4)));
        else
            XYZ = Filament.Data{idx}(:,1:3);
            res(n) = min(pdist2(XYZ,Results(n,3:5)));
        end
    end
end
res = sum(res);

function Ok(hPathsStatsGui)
global Molecule;
global Filament;
global Index;
PathStats=getappdata(hPathsStatsGui.fig,'PathStats');

for n=1:length(PathStats)
    if isreal(Index(n))
        idx = real(Index(n));
        Molecule(idx).PathData = PathStats(n).PathData;
    else
        idx = imag(Index(n));
        Filament(idx).PathData = PathStats(n).PathData;
    end
end
close(hPathsStatsGui.fig);
fShow('Tracks'); %JS Edit to immediately color recolored molecules
    
function Cancel(hPathsStatsGui)
close(hPathsStatsGui.fig);

function Load(hPathsStatsGui)
PathStats=getappdata(hPathsStatsGui.fig,'PathStats');
[FileName, PathName] = uigetfile({'*.mat','FIESTA Path(*.mat)'},'Load FIESTA Path',fShared('GetLoadDir'));
if FileName~=0
    fShared('SetLoadDir',PathName);
    temp_PathStats = fLoad([PathName FileName],'PathStats');
    PathStats = [PathStats temp_PathStats];     
    set(hPathsStatsGui.rLinear,'Value',0);                       
    set(hPathsStatsGui.rAverage,'Value',0);                       
    set(hPathsStatsGui.rPoly2,'Value',0);    
    set(hPathsStatsGui.rPoly3,'Value',0);  
    set(hPathsStatsGui.eAverage,'Enable','off'); 
    set(hPathsStatsGui.tNM,'Enable','off');
    if PathStats(1).AverageDis == -1
        set(hPathsStatsGui.rLinear,'Value',1);                       
    elseif PathStats(1).AverageDis == -2
        set(hPathsStatsGui.rPoly2,'Value',1);   
    elseif PathStats(1).AverageDis == -3
        set(hPathsStatsGui.rPoly3,'Value',1); 
    elseif PathStats(1).AverageDis > 0
        set(hPathsStatsGui.rAverage,'Value',1);          
        set(hPathsStatsGui.eAverage,'Enable','on','String',num2str(PathStats(1).AverageDis)); 
        set(hPathsStatsGui.tNM,'Enable','on');
    end
    set(hPathsStatsGui.lAll,'String',{PathStats.Name},'Value',1,'Max',length(PathStats));  
    set(hPathsStatsGui.cDrift,'Value',PathStats(1).Drift);   
    setappdata(hPathsStatsGui.fig,'PathStats',PathStats);
    Draw(hPathsStatsGui);
end    


function Save(hPathsStatsGui)
PathStats=getappdata(hPathsStatsGui.fig,'PathStats');
Mode=get(gcbo,'UserData');
if strcmp(Mode,'mat');
    [FileName, PathName] = uiputfile({'*.mat','MAT-files (*.mat)'},'Save FIESTA Path',fShared('GetSaveDir'));
    fShared('SetSaveDir',PathName);
    file = [PathName FileName];
    if isempty(findstr('.mat',file))
        file = [file '.mat'];
    end
    hMainGui = getappdata(0,'hMainGui');
    Config = getappdata(hMainGui.fig,'Config'); %#ok<NASGU>
    save(file,'PathStats','Config');
else
    if strcmp(Mode,'single');
        [FileName, PathName] = uiputfile({'*.txt','Delimeted Text (*.txt)'}, 'Save FIESTA Tracks as...',fShared('GetSaveDir'));
        file = [PathName FileName];
        if FileName~=0
            if isempty(findstr('.txt',file))
                file = [file '.txt'];
            end        
            f = fopen(file,'w');
        end
    else
        PathName=uigetdir(fShared('GetSaveDir'));
    end
    if PathName~=0
        fShared('SetSaveDir',PathName);
        for n = 1:length(PathStats)
            if strcmp(Mode,'multiple')
                file=[PathName filesep PathStats(n).Name '.txt'];
                f = fopen(file,'w');
            end
            fprintf(f,'%s - %s%s\n',PathStats(n).Name,PathStats(n).Directory,PathStats(n).File);
            fprintf(f,'frame\ttime[s]\tx-position[nm]\ty-position[nm]\tdistance(to origin)[nm]\tamplitude[ABU]\tpath x-position[nm]\tpath y-Position[nm]\tdistance(along path)[nm]\tsideways(to path)[nm]\n');
            fprintf(f,num2str([PathStats(n).Results(:,1:5) PathStats(n).Results(:,7) PathStats(n).PathData]));
            fprintf(f,'\n'); 
            if strcmp(Mode,'multiple')
                fclose(f);
            end
        end
        if strcmp(Mode,'single')
            fclose(f);
        end
    end
end

function Update(hPathsStatsGui)
global Filament;
global Index;
n=get(hPathsStatsGui.lAll,'Value');
PathStats=getappdata(hPathsStatsGui.fig,'PathStats');
% JS Edit 2022/05/29
% if ~isfield(hPathsStatsGui.pOptions, 'minFrames')
%     minFrames = 6; blinkTolerance = 4;
% else
%     minFrames = getappdata(hPathsStatsGui.pOptions,'minFrames'); blinkTolerance = getappdata(hPathsStatsGui.pOptions,'blinkTolerance');
% end
% End of JS Edit 2022/05/29

if gcbo==hPathsStatsGui.rFilament && ~isreal(Index(n))
    set(hPathsStatsGui.rFilament,'Value',0);  
    return;
end
set(hPathsStatsGui.rLinear,'Value',0);                       
set(hPathsStatsGui.rAverage,'Value',0);                       
set(hPathsStatsGui.rPoly2,'Value',0);    
set(hPathsStatsGui.rPoly3,'Value',0);  
set(hPathsStatsGui.rFilament,'Value',0);  
set(hPathsStatsGui.eAverage,'Enable','off'); 
set(hPathsStatsGui.tRegion,'Enable','off');
set(hPathsStatsGui.tNM,'Enable','off');
if gcbo==hPathsStatsGui.rLinear
    set(hPathsStatsGui.rLinear,'Value',1);                       
elseif gcbo==hPathsStatsGui.rPoly2
    set(hPathsStatsGui.rPoly2,'Value',1);   
elseif gcbo==hPathsStatsGui.rPoly3
    set(hPathsStatsGui.rPoly3,'Value',1); 
elseif gcbo==hPathsStatsGui.rFilament
    set(hPathsStatsGui.rFilament,'Value',1);     
elseif gcbo==hPathsStatsGui.rAverage
    set(hPathsStatsGui.rAverage,'Value',1);          
    set(hPathsStatsGui.eAverage,'Enable','on','String',''); 
    set(hPathsStatsGui.tNM,'Enable','on');
    set(hPathsStatsGui.tRegion,'Enable','on');
    return
end
if gcbo==hPathsStatsGui.lAll||gcbo==hPathsStatsGui.bReset||gcbo==hPathsStatsGui.bDisregard
    i=n(1);
    if gcbo==hPathsStatsGui.bDisregard
        PathStats(i)=[];
        Index(i)=[];
        if isempty(PathStats)
            delete(hPathsStatsGui.fig);
            return;
        end
        if i>length(PathStats)
            i=length(PathStats);
            set(hPathsStatsGui.lAll,'Value',i);
        end
        set(hPathsStatsGui.lAll,'String',{PathStats.Name});
        setappdata(hPathsStatsGui.fig,'PathStats',PathStats);
    end
    if PathStats(i).AverageDis==-1
        set(hPathsStatsGui.rLinear,'Value',1);    
    elseif PathStats(i).AverageDis==-2
        set(hPathsStatsGui.rPoly2,'Value',1);   
    elseif PathStats(i).AverageDis==-3
        set(hPathsStatsGui.rPoly3,'Value',1);
    elseif PathStats(i).AverageDis==-4
        set(hPathsStatsGui.rFilament,'Value',1);        
    elseif PathStats(i).AverageDis>0
        set(hPathsStatsGui.rAverage,'Value',1);          
        set(hPathsStatsGui.eAverage,'Enable','on','String',num2str(PathStats(i).AverageDis)); 
        set(hPathsStatsGui.tNM,'Enable','on');   
    end    
    hPathsStatsGui.Zoom.currentXY = [];
    Draw(hPathsStatsGui);
% JS Edit 2022/06/03 for ranking steppers
elseif gcbo == hPathsStatsGui.bStepSort
    m = str2double(fInputDlg('Expected average stepping rate (nm/frame):', '2.0'));
    neworder = fStepSort(PathStats, m);
    % update list
    tempPathStats = PathStats;
    tempIndex = Index;
    for i = 1:length(PathStats)
        tempPathStats(i) = PathStats(neworder(i));
        tempIndex(i) = Index(neworder(i));
    end
    PathStats = tempPathStats;
    Index = tempIndex;
    
    % Update fig according to the ranking determined in StepSort 
    set(hPathsStatsGui.lAll,'Value',1); % set it to the first
    set(hPathsStatsGui.lAll,'String',{PathStats.Name});
    setappdata(hPathsStatsGui.fig,'PathStats',PathStats);
    
    hPathsStatsGui.Zoom.currentXY = [];
    Draw(hPathsStatsGui);    
% end of JS Edit 2022/06/03
else
    cla(hPathsStatsGui.aPlotXYZ);
    set(hPathsStatsGui.aPlotXYZ,'Visible','off');   
    set(hPathsStatsGui.tCalcPath,'Visible','on');
    drawnow;
    for i = n
        if gcbo==hPathsStatsGui.bAuto
            [param1,resnorm1] = PathFitLinear(PathStats(i).Results(:,3:5));
            [param2,resnorm2] = PathFitCurved(PathStats(i).Results(:,3:5),3); 
            [param3,resnorm3] = PathFitCurved(PathStats(i).Results(:,3:5),4);
            if all(resnorm1*1.1<[resnorm2 resnorm3])
                Path = EvalLinearPath(param1,PathStats(i).Results(:,3:5));
                if i==n(1)
                   set(hPathsStatsGui.rLinear,'Value',1);  
                end
                PathStats(i).AverageDis = -1;
            elseif resnorm2*1.1<resnorm3
                Path = EvalCurvedPath(param2,PathStats(i).Results(:,3:5));
                PathStats(i).AverageDis = -2;
                if i==n(1)
                    set(hPathsStatsGui.rPoly2,'Value',1); 
                end   
            else
                Path = EvalCurvedPath(param3,PathStats(i).Results(:,3:5));
                PathStats(i).AverageDis = -3;
                if i==n(1)
                    set(hPathsStatsGui.rPoly3,'Value',1); 
                end  
            end
        else
            if get(hPathsStatsGui.rLinear,'Value')==1
                [param1,~] = PathFitLinear(PathStats(i).Results(:,3:5));
                Path = EvalLinearPath(param1,PathStats(i).Results(:,3:5));
                PathStats(i).AverageDis=-1;
            elseif get(hPathsStatsGui.rPoly2,'Value')==1
                 [param2,resnorm] = PathFitCurved(PathStats(i).Results(:,3:5),3); 
                 if resnorm<1e+100
                     Path = EvalCurvedPath(param2,PathStats(i).Results(:,3:5));
                     PathStats(i).AverageDis = -2;
                 else
                    Path = EvalLinearPath(param1,PathStats(i).Results(:,3:5));
                    set(hPathsStatsGui.rLinear,'Value',1);  
                    set(hPathsStatsGui.rPoly2,'Value',0);  
                    PathStats(i).AverageDis=-1;
                 end
            elseif get(hPathsStatsGui.rPoly3,'Value')==1
                [param3,resnorm] = PathFitCurved(PathStats(i).Results(:,3:5),4); 
                if resnorm<1e+100
                    Path = EvalCurvedPath(param3,PathStats(i).Results(:,3:5));
                    PathStats(i).AverageDis = -3;
                else
                    Path = EvalLinearPath(param1,PathStats(i).Results(:,3:5));
                    set(hPathsStatsGui.rLinear,'Value',1);  
                    set(hPathsStatsGui.rPoly3,'Value',0);  
                    PathStats(i).AverageDis=-1;
                end
             elseif get(hPathsStatsGui.rFilament,'Value')==1
                 PosFil = getappdata(hPathsStatsGui.fig,'PosFil');      
                 PathFilSelect = getappdata(hPathsStatsGui.fig,'PathFilSelect');      
                 if isempty(PosFil)
                    FilSelect = [Filament.Selected];
                    if all(FilSelect==0)
                        PathFilSelect = ones(size(Filament));
                    else
                        PathFilSelect = zeros(size(Filament));
                        for n=1:length(Filament)
                            if size(Filament(n).Results,1)==1
                                PathFilSelect(n)=1;
                            end
                        end
                    end
                    PathFilSelect = find(PathFilSelect==1);
                    PosFil = InterpolFil(Filament(PathFilSelect));
                    setappdata(hPathsStatsGui.fig,'PosFil',PosFil);
                    setappdata(hPathsStatsGui.fig,'PathFilSelect',PathFilSelect);
                 end
                 nFil = size(PosFil,1);
                 if nFil > 1
                     s = ones(nFil,1)*Inf;
                     for m = 1:nFil
                         s(m) = GetFilament(double(PathStats(i).Results(:,1:5)),Filament(PathFilSelect(m)));
                     end
                     [~,k] = min(s);
                 else
                     k = 1;
                 end
                 Path = EvalFilamentPath(PathStats(i).Results,PosFil(k,:),Filament(PathFilSelect(k)));
                 PathStats(i).AverageDis = -4;
            else
                if ~isempty(get(hPathsStatsGui.eAverage,'String'))
                    AverageDis = round(str2double(get(hPathsStatsGui.eAverage,'String')));
                    Path = AveragePath(PathStats(n).Results(:,1:5),AverageDis);
                    PathStats(i).AverageDis = AverageDis;
                    set(hPathsStatsGui.eAverage,'Enable','on');     
                    set(hPathsStatsGui.rAverage,'Value',1);   
                else
                    Path = [];
                end
            end
        end
        PathStats(i).PathData = Path;
    end

setappdata(hPathsStatsGui.fig,'PathStats',PathStats);
% % JS Edit 2022/05/29
% setappdata(hPathsStatsGui.pOptions,'minFrames',minFrames);
% setappdata(hPathsStatsGui.pOptions,'blinkTolerance',blinkTolerance);
% % End of JS Edit 2022/05/29
Draw(hPathsStatsGui);    

end

% JS Edit 2023/01/09 Crop option
function Crop(hPathsStatsGui)
global Config;
global Molecule;

n = get(hPathsStatsGui.lAll,'Value');
PathStats=getappdata(hPathsStatsGui.fig,'PathStats');
% That is stored in PathStats(n).PathData
xy = PathStats(n).PathData(:,4:5);

% put in NaN for the blank spots. Shouldn't affect the fitting or any post
% processing
frames = PathStats(n).Results(:,1);
frame1 = frames(1); framef = frames(end);
xynew = nan(framef-frame1+1,2);
for i = 1:size(xy,1)
    xynew(frames(i)-frame1+1,:) = xy(i,:);
end

% JS Edit 2023/01/06
% crop option 

[startpoint,endpoint] = parse_pass_v3(xynew);
if endpoint <= startpoint
    error("Select endpoint then startpoint")
end

% only if cropping was actually done should we go through all this madness
if startpoint ~= 1 || endpoint ~= length(xynew)
    
    xynew = xynew(startpoint:endpoint,:);
    % also we will set the new molecule data here so that it is representative
    % on the figure
    
    % Crop your molecule on main fig ... in a very roundabout way to refind it
    s = get(hPathsStatsGui.lAll,'String');
    molidx = [];
    for j = 1:length(Molecule)
        if strcmp(Molecule(j).Name, s{n})
            molidx = [molidx, j];
        end
    end
    % if there is more than one molecule labeled the same, then get the one
    % that is closest
    if length(molidx) > 1
        meandist = 1e10; minj = 0;
        for j = 1:length(molidx)
            mm = molidx(j);
            if size(Molecule(mm).Results,1) == size(PathStats(n).PathData,1)
                dist = mean( sqrt( (Molecule(mm).Results(:,3) - PathStats(n).PathData(:,1)).^2 + (Molecule(mm).Results(:,3) - PathStats(n).PathData(:,2)).^2) );
                if dist < meandist
                    meandist = dist; minj = molidx(j);
                end
            end
        end
    else
        minj = molidx;
    end

    framefill = frame1:framef;
    bidx = find(Molecule(minj).Results(:,1) == framefill(startpoint));
    eidx = find(Molecule(minj).Results(:,1) == framefill(endpoint));

    Molecule(minj).Results = Molecule(minj).Results(bidx:eidx,:);

    % Then update PathStats and truncate accordingly
    PathStats(n).PathData = PathStats(n).PathData(bidx:eidx,:);
    PathStats(n).Results = PathStats(n).Results(bidx:eidx,:);
    if length(PathStats(n).TrackingResults) > length(bidx:eidx)
        PathStats(n).TrackingResults = PathStats(n).TrackingResults(bidx:eidx);
    end
    setappdata(hPathsStatsGui.fig,'PathStats',PathStats);
    Draw(hPathsStatsGui); %then redraw
    
    % Also delete any previous save data on this molecule since it has now
    % changed
    FFolderName = fullfile(Config.Directory{1}, Config.StackName{1}(1:end-6));
    fname = fullfile(FFolderName,strcat(PathStats(n).Name(10:end),'_fiona.mat'));
    if isfile(fname)
        delete([fname])
    end
end



% JS Edit 2022/05/29 to make FIONA in 
function FIONA(hPathsStatsGui)
global Config;
global Molecule;
n = get(hPathsStatsGui.lAll,'Value');
PathStats=getappdata(hPathsStatsGui.fig,'PathStats');
% That is stored in PathStats(n).PathData
xy = PathStats(n).PathData(:,4:5);

% JS Edit 2022/10/05
% Color your molecule ... in a very roundabout way to refind it
s = get(hPathsStatsGui.lAll,'String');
molidx = [];
for j = 1:length(Molecule)
    if strcmp(Molecule(j).Name, s{n})
        molidx = [molidx, j];
    end
end
% if there is more than one molecule labeled the same, then get the one
% that is closest
if length(molidx) > 1
    meandist = 1e10; minj = 0;
    for j = 1:length(molidx)
        mm = molidx(j);
        if size(Molecule(mm).Results,1) == size(PathStats(n).PathData,1)
            dist = mean( sqrt( (Molecule(mm).Results(:,3) - PathStats(n).PathData(:,1)).^2 + (Molecule(mm).Results(:,3) - PathStats(n).PathData(:,2)).^2) );
            if dist < meandist
                meandist = dist; minj = molidx(j);
            end
        end
    end
else
    minj = molidx;
end

Molecule(minj).Color = [0 1 1];

% to save fiona like traces for fiona analysis
% making a separate file between pauses (option)
% with a certain tolerance of blinking

% Below this does separate files for any missing frame, a pause
% frames = PathStats(n).Results(:,1);
% frame1 = frames(1); laststart = 1; ii=1; currDarkFrames = 0; thealphabet = char(97:122);
% 
% minFrames = str2double(fInputDlg('Minimum number of frames to save file:',string(getappdata(hPathsStatsGui.pOptions,'minFrames'))));  
% blinkTolerance = str2double(fInputDlg('Maximum number of dark frames (blink) to concatenate into same file:',string(getappdata(hPathsStatsGui.pOptions,'blinkTolerance')))); 
% setappdata(hPathsStatsGui.pOptions,'minFrames',minFrames); setappdata(hPathsStatsGui.pOptions,'blinkTolerance',blinkTolerance);
% 
% for i = 1:size(xy,1)
%     curr = frame1 + i - laststart + currDarkFrames;
%     if frames(i) ~= curr
%         if frames(i) - curr < blinkTolerance % no abs since should only be positive
%             currDarkFrames = currDarkFrames + frames(i) - curr;
%         else %if not equivalent and not within blinkTolerance, then save as separate file
%             if (i - laststart) >= minFrames
%                 writematrix(xy(laststart:i-1,:)-min(xy(laststart:i-1,:)), strcat(Config.Directory{1},PathStats(n).Name(10:end),thealphabet(ii),'_fiona.txt'), 'Delimiter', 'tab');
%                 ii = ii + 1;
%             end
%             laststart = i;
%             frame1 = frames(i);
%             currDarkFrames = 0;
%         end
%     end
% end
% % make sure to save the end
% if (i - laststart) >= minFrames
%     writematrix(xy(laststart:i,:)-min(xy(laststart:i,:)), strcat(Config.Directory{1},PathStats(n).Name(10:end),thealphabet(ii),'_fiona.txt'), 'Delimiter', 'tab');

% Store in a separate folder for safekeeping
Config.Directory
Config.StackName
FFolderName = fullfile(Config.Directory{1}, Config.StackName{1}(1:end-6));
if ~isfolder(FFolderName)
    mkdir(FFolderName);
end
fname = fullfile(FFolderName,PathStats(n).Name(10:end));

% JS Edit 2022/05/29 to make FIONA
% put in NaN for the blank spots. Shouldn't affect the fitting or any post
% processing
frames = PathStats(n).Results(:,1);
frame1 = frames(1); framef = frames(end);
xynew = nan(framef-frame1+1,2);
for i = 1:size(xy,1)
    xynew(frames(i)-frame1+1,:) = xy(i,:);
end

plotNeighbors = get(hPathsStatsGui.cNeighbor,'Value');
if plotNeighbors == 1
    % JS Edit 2023/01/31
    % we should extend the path and more accurately plot the neighbors onto the
    % x,y coordinates of the channel we have determined
    
    % first we need to find the least squares path if it isn't empty, which
    % usually is cubic which may cause problems
    if (PathStats(n).AverageDis == 0) || (PathStats(n).AverageDis == -5) 
        % use path data for the saved version
        [~,resnorm1] = PathFitLinear(PathStats(n).PathData(:,1:3));
        [~,resnorm2] = PathFitCurved(PathStats(n).PathData(:,1:3),3);
        [~,resnorm3] = PathFitCurved(PathStats(n).PathData(:,1:3),4);
        if all(resnorm1*1.1<[resnorm2 resnorm3])
            PathStats(n).AverageDis = -1;
        elseif resnorm2*1.1<resnorm3
            PathStats(n).AverageDis = -2;
        else
            PathStats(n).AverageDis = -3;
        end
    end
    
    % now we need to do something to extend the path to make it easier for
    % the neighbors to find
    % first, make an interpolation on latent variable so that it approximately goes a nm (we'll stick with 2d for now)
    path_extend = 200; %nm

    % interpolate between each data point by 1 nm
    Delta = PathStats(n).PathData(2:end,1:4)-PathStats(n).PathData(1:end-1,1:4);
    InterpResults = PathStats(n).PathData(1,1:3); % initialize with first point
    for i = 1:size(Delta,1)
        dstep = round(Delta(i,4)); % approximately 1 nm per step
        step = Delta(i,1:3)/dstep;
        xi = PathStats(n).PathData(i,1)*ones(1,dstep) + step(1)*(1:dstep);
        yi = PathStats(n).PathData(i,2)*ones(1,dstep) + step(2)*(1:dstep);
        zi = PathStats(n).PathData(i,3)*ones(1,dstep) + step(3)*(1:dstep);
        InterpResults = [InterpResults; xi' yi' zi'];
    end
    % Now we want to do the ends, a linear extrapolation at the end
    % this isn't as good for high curves but deals with edge cases
    w = 1:min(size(InterpResults,1),path_extend);
    % beginning fits
    xbeginmdl = fit(w',InterpResults(1:length(w),1),'linear');
    ybeginmdl = fit(w',InterpResults(1:length(w),2),'linear');
    % ending fits
    
    xendmdl = fit(w',InterpResults(end-length(w)+1:end,1),'linear');
    yendmdl = fit(w',InterpResults(end-length(w)+1:end,2),'linear');
    if all(~isnan(InterpResults(:,3)))
        zbeginmdl = fit(w',InterpResults(1:length(w),3),'linear');
        zendmdl = fit(w',InterpResults(end-length(w):end,3),'linear');
        zb = zbeginmdl(-path_extend+1:0)';
        ze = zendmdl(length(w)+1:length(w)+path_extend)';
    else
        zb = nan(length(w),1);
        ze = nan(length(w),1);
    end
    tb = -path_extend+1:0;
    te = length(w)+1:length(w)+path_extend;
    
    % Now fill into InterpResults in proper places
    InterpResults = [xbeginmdl(tb), ybeginmdl(tb), zb; InterpResults];
    InterpResults = [InterpResults; xendmdl(te), yendmdl(te), ze];
    
    
    % then, make a path with interpolation to place the neighbors
    if PathStats(n).AverageDis == -1
        [param1,~] = PathFitLinear(PathStats(n).Results(:,3:5));
        InterpPath = EvalLinearPath(param1,InterpResults);
    elseif PathStats(n).AverageDis == -2
        [param2,~] = PathFitCurved(PathStats(n).Results(:,3:5),3);
        InterpPath = EvalCurvedPath(param2,InterpResults);
    elseif PathStats(n).AverageDis == -3
        [param3,~] = PathFitCurved(PathStats(n).Results(:,3:5),4);
        InterpPath = EvalCurvedPath(param3,InterpResults);
    else % can't extend path, fit is individual
        InterpPath = PathStats(n).PathData;
    end
    
    neighbors = findNeighbors(InterpPath(:,1:2));
    neighbor_txy = cell(length(neighbors),1);
    neighbor_exist_thresh = 10;
    
    for m = length(neighbors):-1:1
        Res = Molecule(neighbors(m)).Results;
        if size(Res,1) < neighbor_exist_thresh
            neighbor_txy(m) = [];
        else
        txy_reorient = nan( Res(end,1)-Res(1,1) ,3);
        B = [PathStats(n).Results(end,3) - PathStats(n).Results(1,3), PathStats(n).Results(end,4) - PathStats(n).Results(1,4),0]; %for cross product
        for p = 1:size(Res,1)
%             % Plot on direct path of molecule
%             [npos, nidx] = min( sqrt( (PathStats(n).PathData(:,1) - Res(p,3)).^2 + (PathStats(n).PathData(:,2) - Res(p,4)).^2) );
%             mpos = norm (PathStats(n).PathData(nidx,1:2) - PathStats(n).PathData(1,1:2));
            % Plot on interpolated path
            [npos, nidx] = min( sqrt( (InterpPath(:,1) - Res(p,3)).^2 + (InterpPath(:,2) - Res(p,4)).^2) );
            mpos = norm (InterpPath(nidx,1:2) - PathStats(n).PathData(1,1:2));
            
            % to get correct sign off axis, use a cross product on vectors
            % for overall direction
            C = [Res(p,3) - InterpPath(nidx,1), Res(p,4) - InterpPath(nidx,2),0]; % to nearest point
            D = [Res(p,3) - PathStats(n).Results(1,3), Res(p,4) - PathStats(n).Results(1,4),0]; % to beginning
            BXC = cross(B,C); npos = npos * sign(BXC(3));
            theta = acos(dot(B,D)/norm(B)/norm(D)); mpos = mpos * sign(cos(theta));
            
            txy_reorient(Res(p,1)-Res(1,1)+1,:) = [Res(p,1) - frame1 + 1, mpos, npos];
        end
        neighbor_txy{m,1} = txy_reorient;
        end
    end
else
    neighbor_txy = {};
end

% New way to save
if ~isfile(strcat(fname,'_fiona.mat')) % Don't overwrite existing file
data.xy = xynew;
data.yx = xynew(:,[2,1]);
data.neighbors = neighbor_txy;
data.time = PathStats(n).Results(:,2);
save(strcat(fname,'_fiona.mat'),'data');
end

% at the end run FIONAviewer
FIONAviewer(fullfile(FFolderName,strcat(PathStats(n).Name(10:end),'_fiona.mat')))


% JS Edit 2022/06/03 to compile plots with the click of a button
function StepStats(hPathsStatsGui)
fStepStats()


function Draw(hPathsStatsGui)
cline = [0.4, 0.5, 1];
set(hPathsStatsGui.tCalcPath,'Visible','off');
set(hPathsStatsGui.aPlotXYZ,'Visible','on');
set(hPathsStatsGui.aPlotXYZ,'Color',get(hPathsStatsGui.fig,'Color'));
PathStats=getappdata(hPathsStatsGui.fig,'PathStats');
idx=get(hPathsStatsGui.lAll,'Value');
idx=idx(1);
XPlotPath = PathStats(idx).PathData(:,1)-min(PathStats(idx).Results(:,3));
YPlotPath = PathStats(idx).PathData(:,2)-min(PathStats(idx).Results(:,4));
Dis = real(PathStats(idx).PathData(:,4));
XPlot = PathStats(idx).Results(:,3)-min(PathStats(idx).Results(:,3));
YPlot = PathStats(idx).Results(:,4)-min(PathStats(idx).Results(:,4));
if any(isnan(PathStats(idx).PathData(:,3)))
    ZPlot = [];
    ZPlotPath = [];
else
    ZPlotPath = PathStats(idx).PathData(:,3)-min(PathStats(idx).Results(:,5));
    ZPlot = PathStats(idx).Results(:,5)-min(PathStats(idx).Results(:,5));
end
if ((max(XPlot)-min(XPlot)) > 5000) || ((max(YPlot)-min(YPlot)) > 5000)
    scale=1000;
    units=['[' char(181) 'm]'];
else
    scale=1;
    units='[nm]';
end
data = [Dis XPlotPath YPlotPath ZPlotPath];
data = sortrows(data,1);
XPlotPath = data(:,2);
YPlotPath = data(:,3);
if isempty(ZPlot)
    plot(hPathsStatsGui.aPlotXYZ,XPlot/scale,YPlot/scale,'Color',cline,'LineStyle','-','Marker','*');
    line(hPathsStatsGui.aPlotXYZ,XPlotPath/scale,YPlotPath/scale,'Color','green','LineStyle','-','Marker','none');
    zoom(hPathsStatsGui.aPlotXYZ,'on');
    axis equal
else
    units='[nm]';
    ZPlotPath = data(:,4);
    scatter3(hPathsStatsGui.aPlotXYZ,XPlot,YPlot,ZPlot,'MarkerEdgeColor','k','MarkerFaceColor',[0 .75 .75]);
    hold(hPathsStatsGui.aPlotXYZ,'on');
    plot3(hPathsStatsGui.aPlotXYZ,XPlotPath,YPlotPath,ZPlotPath,'Color','red','LineStyle','-','Marker','none','LineWidth',5);
    hold(hPathsStatsGui.aPlotXYZ,'off');
    rotate3d(hPathsStatsGui.aPlotXYZ,'on');
end
set(hPathsStatsGui.aPlotXYZ,'YDir','reverse');
xlabel(hPathsStatsGui.aPlotXYZ,['X-Position ' units],'Color','White');
ylabel(hPathsStatsGui.aPlotXYZ,['Y-Position ' units],'COlor','White');
set(hPathsStatsGui.aPlotXYZ,'XColor','white'); set(hPathsStatsGui.aPlotXYZ,'YColor','white');
if ~isempty(ZPlot)
    zlabel(hPathsStatsGui.aPlotXYZ,'z-Position [nm]','Color','white');
end
set(hPathsStatsGui.aPlotXYZ,'Color',get(hPathsStatsGui.fig,'Color'));
XPlot = PathStats(idx).Results(:,2);
YPlot = (Dis-Dis(1));
YPlotOld = PathStats(idx).Results(:,6);
if (max(YPlot)-min(YPlot)) > 5000 || (max(YPlotOld)-min(YPlotOld)) > 5000
    yscale=1000;
    units=['[' char(181) 'm]'];
else
    yscale=1;
    units='[nm]';
end 
% JS Edit 2023/01/09 for crop plotting
%plot(hPathsStatsGui.aPlotDist,XPlot,YPlotOld/yscale,'Color','magenta','LineStyle','--','Marker','none');
plot(hPathsStatsGui.aPlotDist,XPlot,(YPlotOld-YPlotOld(1))/yscale,'Color','magenta','LineStyle','--','Marker','none');
line(hPathsStatsGui.aPlotDist,XPlot,YPlot/yscale,'Color',cline,'LineStyle','-','Marker','none');
xlabel(hPathsStatsGui.aPlotDist,'Time [sec]','Color','white');
ylabel(hPathsStatsGui.aPlotDist,['Distance along path ' units],'Color','white');
set(hPathsStatsGui.aPlotDist,'XColor','white'); set(hPathsStatsGui.aPlotDist,'YColor','white'); 
if ~any(isnan(PathStats(idx).PathData(:,6)))
    yyaxis left;
end
set(hPathsStatsGui.aPlotDist,'Color',get(hPathsStatsGui.fig,'Color'));
XPlot = PathStats(idx).Results(:,2);
YPlot = PathStats(idx).PathData(:,5);
plot(hPathsStatsGui.aPlotSide,XPlot,YPlot,'Color',cline,'LineStyle','-','Marker','none');
xlabel(hPathsStatsGui.aPlotSide,'Time [sec]','Color','white');
ylabel(hPathsStatsGui.aPlotSide,'Sideways motion [nm]','Color','white');
set(hPathsStatsGui.aPlotSide,'XColor','white'); set(hPathsStatsGui.aPlotSide,'YColor','white'); 
if ~any(isnan(PathStats(idx).PathData(:,6)))
    yyaxis right;
    XPlot = PathStats(idx).Results(:,2);
    YPlot = PathStats(idx).PathData(:,6);
    plot(hPathsStatsGui.aPlotSide,XPlot,YPlot,'Color','red','LineStyle','-','Marker','none');
    ylabel(hPathsStatsGui.aPlotSide,'Height motion [nm]');   
end
set(hPathsStatsGui.aPlotSide,'Color',get(hPathsStatsGui.fig,'Color'));
setappdata(0,'hPathsStatsGui',hPathsStatsGui);

function res=Linear2D(param,p,XY)
b = [p(2,1)-p(1,1) param(2,1)-param(1,1)];
pa = [XY(:,1)-p(1,1) XY(:,2)-param(1,1)];
res = (b(1).*pa(:,2) - b(2).*pa(:,1)) / sqrt( b(1)^2 + b(2)^2);

function res=Linear3D(param,p,XYZ)
b = [p(2,1)-p(1,1) param(2,1)-param(1,1) param(2,2)-param(1,2)];
pa = [XYZ(:,1)-p(1,1) XYZ(:,2)-param(1,1) XYZ(:,3)-param(1,2)];
d = [(b(2).*pa(:,3) - b(3).*pa(:,2)) (b(3).*pa(:,1) - b(1).*pa(:,3)) (b(1).*pa(:,2) - b(2).*pa(:,1))];
res = sqrt(sum((d').^2)) / sqrt( b(1)^2 + b(2)^2 + b(3)^2);

function [param,resnorm]=PathFitLinear(XYZ)
XYZ = double(XYZ);
flip = 0;
if max(XYZ(:,2))-min(XYZ(:,2))>max(XYZ(:,1))-min(XYZ(:,1))
    XYZ(:,1:2)=fliplr(XYZ(:,1:2));
    flip = 1;
end
[param0(1,1),kmin] = min(XYZ(:,1));
[param0(2,1),kmax]=max(XYZ(:,1));
param0(1,2)=min(XYZ(kmin,2));
param0(2,2)=max(XYZ(kmax,2));
if all(~isnan(XYZ(:,3)))
    param0(1,3)=min(XYZ(kmin,3));
    param0(2,3)=max(XYZ(kmax,3));
end
options = optimset('lsqnonlin');
options.Display = 'off';
param = param0(:,2:end);
try
    if size(param0,2)==2
        [param,resnorm] = lsqnonlin(@Linear2D,param,[],[],options,param0(:,1),XYZ(:,1:2)); 
    else
        [param,resnorm] = lsqnonlin(@Linear3D,param,[],[],options,param0(:,1),XYZ); 
    end
    param = [param0(:,1) param];
catch %#ok<CTCH>
    param=param0;
    resnorm=1e100;
end
if flip==1
    param(:,1:2) = fliplr(param(:,1:2));
end
if XYZ(1,1)>XYZ(end,1)
    param = flipud(param);
end

function Path = EvalLinearPath(param,XYZ)
if size(param,2)==2
    param(:,3) = 0;
    XYZ(:,3) = 0;
end
XYZ = double(XYZ);
b = [param(2,1)-param(1,1) param(2,2)-param(1,2) param(2,3)-param(1,3)];
pa = [XYZ(:,1)-param(1,1) XYZ(:,2)-param(1,2) XYZ(:,3)-param(1,3)];
t0 = (pa(:,1)*b(1) + pa(:,2)*b(2) + pa(:,3)*b(3)) / (b(1)^2 + b(2)^2 + b(3)^2);
Path(:,1:3) = [param(1,1)+t0*b(1) param(1,2)+t0*b(2) param(1,3)+t0*b(3)];
dis = sqrt( (t0*b(1)).^2 + (t0*b(2)).^2  + (t0*b(3)).^2).*sign(t0);
dis = dis-dis(1);
if mean(dis)<0
    dis = -dis;
end
Path(:,4) = dis;
norm_b = sqrt(sum(b.^2));
nb = [b(:,1)./norm_b b(:,2)./norm_b b(:,3)./norm_b];
Rz = zeros(3,3);
Rz(1,1) = nb(1);
Rz(2,2) = nb(1);
Rz(1,2) = nb(2);
Rz(2,1) = -nb(2);
Rz(3,3) = 1;
nnb = (Rz*nb')';
a = XYZ-Path(:,1:3);
na = zeros(size(a));
na(:,1) = Rz(1,1)*a(:,1) + Rz(1,2)*a(:,2) + Rz(1,3)*a(:,3);
na(:,2) = Rz(2,1)*a(:,1) + Rz(2,2)*a(:,2) + Rz(2,3)*a(:,3);
na(:,3) = Rz(3,1)*a(:,1) + Rz(3,2)*a(:,2) + Rz(3,3)*a(:,3);
Ry = zeros(3,3);
Ry(1,1) = nnb(1);
Ry(3,3) = nnb(1);
Ry(1,3) = nnb(3);
Ry(3,1) = -nnb(3);
Ry(2,2) = 1;
X(:,1) = Ry(1,1)*na(:,1) + Ry(1,2)*na(:,2) + Ry(1,3)*na(:,3);
X(:,2) = -( Ry(2,1)*na(:,1) + Ry(2,2)*na(:,2) + Ry(2,3)*na(:,3));
X(:,3) = Ry(3,1)*na(:,1) + Ry(3,2)*na(:,2) + Ry(3,3)*na(:,3);
Path(:,5:6) = X(:,2:3);
if all(Path(:,3)==0)
    Path(:,3) = NaN;
    Path(:,6) = NaN;
end

function res = getCurvedPath2D(Path,XY)
seglen = sqrt(sum(diff(Path,[],1).^2,2));
t0 = [0;cumsum(seglen)/sum(seglen)];
spl = csapi(t0,Path');
tt = -0.2:1/(sum(seglen)*1.4):1.2;
V = fnval(spl,tt)';
%res = zeros(1,N);
res = min(pdist2(V,XY));

function res = getCurvedPath3D(Path,XYZ)
%N = size(XYZ,1);
seglen = sqrt(sum(diff(Path,[],1).^2,2));
t0 = [0;cumsum(seglen)/sum(seglen)];
spl = csapi(t0,Path');
tt = -0.2:1/(sum(seglen)*1.4*0.2):1.2;
V = fnval(spl,tt)';
%res = zeros(1,N);
res = min(pdist2(V,XYZ));

function res=Curved2D(params,Path,XY)
n = size(Path,1);
Path(2:end-1,1) = params(1:n-2);
Path(:,2) = params(n-1:2*n-2);
if Path(2,1)>Path(3,1)
    Path(2:3,:) = flipud(Path(2:3,:));
end
res = getCurvedPath2D(Path,XY);

function res=Curved3D(params,Path,XYZ)
n = size(Path,1);
Path(2:end-1,1) = params(1:n-2);
Path(:,2) = params(n-1:2*n-2);
Path(:,3) = params(2*n-1:3*n-2);
if Path(2,1)>Path(3,1)
    Path(2:3,:) = flipud(Path(2:3,:));
end
res = getCurvedPath3D(Path,XYZ);

function [path,resnorm]=PathFitCurved(XYZ,knots)
XYZ = double(XYZ);
flip = 0;
if max(XYZ(:,2))-min(XYZ(:,2))>max(XYZ(:,1))-min(XYZ(:,1))
    XYZ(:,1:2)=fliplr(XYZ(:,1:2));
    flip = 1;
end
[param0(1,1),kmin] = min(XYZ(:,1));
[param0(2,1),kmax]=max(XYZ(:,1));
param0(1,2)=min(XYZ(kmin,2));
param0(2,2)=max(XYZ(kmax,2));
if all(~isnan(XYZ(:,3)))
    param0(1,3)=min(XYZ(kmin,3));
    param0(2,3)=max(XYZ(kmax,3));
end
param0(knots,:) = param0(end,:);
for n = 2:knots-1
    param0(n,:) = ((knots-n)*param0(1,:)+(n-1)*param0(knots,:))/(knots-1);
    [~,k] = min(abs(XYZ(:,1)-param0(n,1)));
    param0(n,2) = XYZ(k,2);
    if ~isnan(XYZ(k,3))
        param0(n,3) = XYZ(k,3);    
    end
end
lb = param0;
ub = param0;
lb(2:end-1,1) = min(lb(:,1));
ub(2:end-1,1) = max(ub(:,1));
lb(:,2:end) = -Inf;
ub(:,2:end) = +Inf;
param = param0;
param(1)=[];
param(knots-1)=[];
lb(1)=[];
lb(knots-1)=[];
ub(1)=[];
ub(knots-1)=[];
options = optimset('lsqnonlin');
options.Display = 'off';
path = param0;
try
    if size(param0,2)==2
        [param,resnorm]= lsqnonlin(@Curved2D,param,lb,ub,options,param0,XYZ(:,1:2));
        path(2:end-1,1) = param(1:knots-2);
        path(:,2) = param(knots-1:2*knots-2);
    else
        [param,resnorm]= lsqnonlin(@Curved3D,param,lb,ub,options,param0,XYZ);
        path(2:end-1,1) = param(1:knots-2);
        path(:,2) = param(knots-1:2*knots-2);
        path(:,3) = param(2*knots-1:3*knots-2);
    end
catch %#ok<CTCH>
    path=param0;
    resnorm=1e100;
end
if flip==1
    path(:,1:2) = fliplr(path(:,1:2));
end
if XYZ(1,1)>XYZ(end,1)
    path = flipud(path);
end

function Path = EvalFilamentPath(Results,PosFil,Filament)
N = size(Results,1);
idx = zeros(N,1);
XYZ = Results(:,3:5);
V = zeros(N,3);
b = zeros(N,3);
for n = 1:N
    if length(PosFil)>1
        [~,k] = min(abs(Results(n,1)-Filament.Results(:,1)));
        frame = k(1);
        PXYZ = PosFil{frame};
        if size(PXYZ,2)==2
            PXYZ(:,3) = 0;
            XYZ(:,3) = 0;
        end
        R = sqrt((PXYZ(:,1)-XYZ(n,1)).^2+(PXYZ(:,2)-XYZ(n,2)).^2+(PXYZ(:,3)-XYZ(n,3)).^2);
        [~,m] = min(R);
        V(n,1:3) = PXYZ(m,1:3);
        if m>=length(R)
            m=length(R)-1;
        end
        b(n,1:3) = PXYZ(m+1,1:3)-PXYZ(m,1:3);
        idx(n) = n;
    else
        V = PosFil{1};
        if size(V,2)==2
            V(:,3) = 0;
            XYZ(:,3) = 0;
        end
        b = zeros(size(V));
        R = sqrt((V(:,1)-XYZ(n,1)).^2+(V(:,2)-XYZ(n,2)).^2+(V(:,3)-XYZ(n,3)).^2);
        [~,idx(n)] = min(R);
    end
end
  
segpath = [0;cumsum(sqrt(sum(diff(V,[],1).^2,2)))];
dis = segpath(idx);
if length(PosFil)==1
    dis = dis-segpath(5001) + (segpath(end-5000)-dis)*1i;
end

Path(:,1:3) = V(idx,:);
Path(:,4) = dis;

a = XYZ - V(idx,:);
if all(all(b==0))
    b(2:end,:) = V(2:end,:) - V(1:end-1,:);
end

%get norm vector in direction of path
norm_b = sqrt(sum((b').^2))';
nb = [b(:,1)./norm_b b(:,2)./norm_b b(:,3)./norm_b];
nb(isnan(nb))=0;
%get rotation matrix for roation around z axis 
Rz = zeros(3,3,size(V,1));
Rz(1,1,:) = nb(:,1);
Rz(2,2,:) = nb(:,1);
Rz(1,2,:) = nb(:,2);
Rz(2,1,:) = -nb(:,2);
Rz(3,3,:) = 1;

nnb = zeros(size(V));
nnb(:,1) = squeeze(Rz(1,1,:)).*nb(:,1) + squeeze(Rz(1,2,:)).*nb(:,2) + squeeze(Rz(1,3,:)).*nb(:,3);
nnb(:,2) = squeeze(Rz(2,1,:)).*nb(:,1) + squeeze(Rz(2,2,:)).*nb(:,2) + squeeze(Rz(2,3,:)).*nb(:,3);
nnb(:,3) = squeeze(Rz(3,1,:)).*nb(:,1) + squeeze(Rz(3,2,:)).*nb(:,2) + squeeze(Rz(3,3,:)).*nb(:,3);

na = zeros(size(a));
na(:,1) = squeeze(Rz(1,1,idx)).*a(:,1) + squeeze(Rz(1,2,idx)).*a(:,2) + squeeze(Rz(1,3,idx)).*a(:,3);
na(:,2) = squeeze(Rz(2,1,idx)).*a(:,1) + squeeze(Rz(2,2,idx)).*a(:,2) + squeeze(Rz(2,3,idx)).*a(:,3);
na(:,3) = squeeze(Rz(3,1,idx)).*a(:,1) + squeeze(Rz(3,2,idx)).*a(:,2) + squeeze(Rz(3,3,idx)).*a(:,3);

Ry = zeros(3,3,size(V,1));
Ry(1,1,:) = nnb(:,1);
Ry(3,3,:) = nnb(:,1);
Ry(1,3,:) = nnb(:,3);
Ry(3,1,:) = -nnb(:,3);
Ry(2,2,:) = 1;

X = zeros(size(na));
X(:,1) = squeeze(Ry(1,1,idx)).*na(:,1) + squeeze(Ry(1,2,idx)).*na(:,2) + squeeze(Ry(1,3,idx)).*na(:,3);
X(:,2) = - (squeeze(Ry(2,1,idx)).*na(:,1) + squeeze(Ry(2,2,idx)).*na(:,2) + squeeze(Ry(2,3,idx)).*na(:,3));
X(:,3) = squeeze(Ry(3,1,idx)).*na(:,1) + squeeze(Ry(3,2,idx)).*na(:,2) + squeeze(Ry(3,3,idx)).*na(:,3);
Path(:,5:6) = X(:,2:3);
if all(Path(:,3)==0)
    Path(:,3) = NaN;
    Path(:,6) = NaN;
end

function Path = EvalCurvedPath(param,XYZ)
if size(param,2)==2
    param(:,3) = 0;
    XYZ(:,3) = 0;
end
XYZ = double(XYZ);
seglen = sqrt(sum(diff(param,[],1).^2,2));
t0 = [0;cumsum(seglen)/sum(seglen)];
spl = csapi(t0,param');
tt = -0.2:1/(sum(seglen)*1.4*10):1.2;
V = fnval(spl,tt)';
N = size(XYZ,1);
idx = zeros(N,1);
for n = 1:N
    R = sqrt((V(:,1)-XYZ(n,1)).^2+(V(:,2)-XYZ(n,2)).^2+(V(:,3)-XYZ(n,3)).^2);
    [~,idx(n)] = min(R);
end
segpath = [0;cumsum(sqrt(sum(diff(V,[],1).^2,2)))];
dis = segpath(idx);
dis = dis-dis(1);
if mean(dis)<0
    dis = -dis;
end

Path(:,1:3) = V(idx,:);
Path(:,4) = dis;

a = XYZ - V(idx,:);

b = zeros(size(V));
b(2:end,:) = V(2:end,:) - V(1:end-1,:);

%get norm vector in direction of path
norm_b = sqrt(sum((b').^2))';
nb = [b(:,1)./norm_b b(:,2)./norm_b b(:,3)./norm_b];
nb(isnan(nb))=0;
%get rotation matrix for roation around z axis 
Rz = zeros(3,3,size(V,1));
Rz(1,1,:) = nb(:,1);
Rz(2,2,:) = nb(:,1);
Rz(1,2,:) = nb(:,2);
Rz(2,1,:) = -nb(:,2);
Rz(3,3,:) = 1;

nnb = zeros(size(V));
nnb(:,1) = squeeze(Rz(1,1,:)).*nb(:,1) + squeeze(Rz(1,2,:)).*nb(:,2) + squeeze(Rz(1,3,:)).*nb(:,3);
nnb(:,2) = squeeze(Rz(2,1,:)).*nb(:,1) + squeeze(Rz(2,2,:)).*nb(:,2) + squeeze(Rz(2,3,:)).*nb(:,3);
nnb(:,3) = squeeze(Rz(3,1,:)).*nb(:,1) + squeeze(Rz(3,2,:)).*nb(:,2) + squeeze(Rz(3,3,:)).*nb(:,3);

na = zeros(size(a));
na(:,1) = squeeze(Rz(1,1,idx)).*a(:,1) + squeeze(Rz(1,2,idx)).*a(:,2) + squeeze(Rz(1,3,idx)).*a(:,3);
na(:,2) = squeeze(Rz(2,1,idx)).*a(:,1) + squeeze(Rz(2,2,idx)).*a(:,2) + squeeze(Rz(2,3,idx)).*a(:,3);
na(:,3) = squeeze(Rz(3,1,idx)).*a(:,1) + squeeze(Rz(3,2,idx)).*a(:,2) + squeeze(Rz(3,3,idx)).*a(:,3);

Ry = zeros(3,3,size(V,1));
Ry(1,1,:) = nnb(:,1);
Ry(3,3,:) = nnb(:,1);
Ry(1,3,:) = nnb(:,3);
Ry(3,1,:) = -nnb(:,3);
Ry(2,2,:) = 1;

X = zeros(size(na));
X(:,1) = squeeze(Ry(1,1,idx)).*na(:,1) + squeeze(Ry(1,2,idx)).*na(:,2) + squeeze(Ry(1,3,idx)).*na(:,3);
X(:,2) = -( squeeze(Ry(2,1,idx)).*na(:,1) + squeeze(Ry(2,2,idx)).*na(:,2) + squeeze(Ry(2,3,idx)).*na(:,3));
X(:,3) = squeeze(Ry(3,1,idx)).*na(:,1) + squeeze(Ry(3,2,idx)).*na(:,2) + squeeze(Ry(3,3,idx)).*na(:,3);
Path(:,5:6) = X(:,2:3);
if all(Path(:,3)==0)
    Path(:,3) = NaN;
    Path(:,6) = NaN;
end

function Path = AveragePath(Results,DisRegion)
nData=size(Results,1);
XYZ= double(Results(:,3:5));
if any(isnan(XYZ(:,3)))
    XYZ(:,3) = 0;    
end
p=2;
n=1;
param = XYZ(1,:);
while n<=nData
    IN = sqrt( ( XYZ(n,1)-XYZ(:,1)).^2 + ( XYZ(n,2)-XYZ(:,2)).^2 + ( XYZ(n,3)-XYZ(:,3)).^2 ) <DisRegion;
    if sum(IN)==1
        if n>1 && n<nData
            param(p,:) = XYZ(n,:);
            p = p+1;
        end
        n = n+1;
    else
        k = find(~IN);
        k_start = k(find(k<n,1,'last'))+1;
        k_end = k(find(k>n,1,'first'))-1;
        if isempty(k_start)
            k_start=1;
        end
        if isempty(k_end)
            k_end=nData;
        end
        param(p,:) = mean(XYZ(k_start:k_end,:));
        n = k_end+1;
        p=p+1;
    end
end
param(end+1,:) = XYZ(end,:);
Path = EvalCurvedPath(param,XYZ);

%JS Edit 2022/09/26 for finding nearby molecules
function nbs = findNeighbors(MolXY, CheckChannel)
global Molecule;

if nargin < 2
    CheckChannel = 2;
end

% Copied from fRightPanel NewScan
nX=MolXY(:,1);
nY=MolXY(:,2);
ScanSize=100; % in nm, make it so that they share a point within a pixel
d=[0; cumsum(sqrt((nX(2:end)-nX(1:end-1)).^2 + (nY(2:end)-nY(1:end-1)).^2))];
dt=max(d)/round(max(d));
id=(0:round(max(d)))'*dt;
scan_length=length(id);
idx = nearestpoint(id,d);
X=zeros(scan_length,1);
Y=zeros(scan_length,1);
dis = id-d(idx);
dis(1)=0;
dis(end)=0;
X(dis==0) = nX(idx(dis==0));
Y(dis==0) = nY(idx(dis==0));
X(dis>0) = nX(idx(dis>0))+(nX(idx(dis>0)+1)-nX(idx(dis>0)))./(d(idx(dis>0)+1)-d(idx(dis>0))).*dis(dis>0);
Y(dis>0) = nY(idx(dis>0))+(nY(idx(dis>0)+1)-nY(idx(dis>0)))./(d(idx(dis>0)+1)-d(idx(dis>0))).*dis(dis>0);
X(dis<0) = nX(idx(dis<0))+(nX(idx(dis<0)-1)-nX(idx(dis<0)))./(d(idx(dis<0)-1)-d(idx(dis<0))).*dis(dis<0);
Y(dis<0) = nY(idx(dis<0))+(nY(idx(dis<0)-1)-nY(idx(dis<0)))./(d(idx(dis<0)-1)-d(idx(dis<0))).*dis(dis<0);
iX=zeros(2*ScanSize+1,scan_length);
iY=zeros(2*ScanSize+1,scan_length);
n=zeros(scan_length,3);
for i=1:length(X)
    if i==1   
        v=[X(i+1)-X(i) Y(i+1)-Y(i) 0];
        n(i,:)=[v(2) -v(1) 0]/norm(v); 
    elseif i==length(X)
        v=[X(i)-X(i-1) Y(i)-Y(i-1) 0];
        n(i,:)=[v(2) -v(1) 0]/norm(v);
    else
        v1=[X(i+1)-X(i) Y(i+1)-Y(i) 0];
        v2=-[X(i)-X(i-1) Y(i)-Y(i-1) 0];
        n(i,:)=v1/norm(v1)+v2/norm(v2); 
        if norm(n(i,:))==0
            n(i,:)=[v1(2) -v1(1) 0]/norm(v1);
        else
            n(i,:)=n(i,:)/norm(n(i,:));
        end
        z=cross(v1,n(i,:));
        if z(3)>0
            n(i,:)=-n(i,:);
        end
    end
    iX(:,i)=linspace(X(i)+ScanSize*n(i,1),X(i)-ScanSize*n(i,1),2*ScanSize+1)';
    iY(:,i)=linspace(Y(i)+ScanSize*n(i,2),Y(i)-ScanSize*n(i,2),2*ScanSize+1)';
end

iX = reshape(iX,numel(iX),1);
iY = reshape(iY,numel(iY),1);
% find boundary for our polygon
bd = boundary(iX,iY);

% check if in polygon for all molecules in other channel (rudimentary, but
% fair)
nbs = [];
for j = 1:length(Molecule)
    if Molecule(j).Channel == CheckChannel %only do if in opposite channel
        inbd = inpolygon(Molecule(j).Results(:,3), Molecule(j).Results(:,4), iX(bd), iY(bd));
        if any(inbd > 0)
            nbs = [nbs, j];
        end
    end
end
