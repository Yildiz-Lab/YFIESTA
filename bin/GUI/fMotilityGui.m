function fMotilityGui(func,varargin)
% JS 2022/07/14 For Yildiz Motility Assays.
switch func
    case 'Create',Create(varargin{:});
    case 'Draw',Draw(varargin{1});
    case 'SaveFig',SaveFig(varargin{1});
    case 'ExportCSV',ExportCSV(varargin{1}); % Preserve the existing callback name.
    case 'Cancel',close(varargin{1}.fig);
end

function Create(varargin)
global Molecule Config TimeInfo Stack
hMainGui=getappdata(0,'hMainGui');
if isempty(Molecule),fMsgDlg('No molecules are loaded!','error');return;end
if isempty(varargin),options=fMotilityStatsDlg(numel(hMainGui.Region));else,options=varargin{1};end
if isempty(options),return;end
if ~options.pool&&~any([Molecule.Selected]==1),fMsgDlg('No molecules selected!','error');return;end
try
    channel=hMainGui.Values.FrameIdx(1);interval=Config.Time(min(channel,numel(Config.Time)));
    context=struct('PixelSizeNm',hMainGui.Values.PixSize,'Channel',channel,'FrameIntervalMs',interval, ...
        'FrameRange',[str2double(get(hMainGui.RightPanel.pTools.eKymoStart,'String')) ...
        str2double(get(hMainGui.RightPanel.pTools.eKymoEnd,'String'))],'TimeMs',[],'Scan',hMainGui.Scan);
    if numel(TimeInfo)>=channel,context.TimeMs=TimeInfo{channel};end
    report=fYildizMotility('analyze',Molecule,hMainGui.Region,context,options);
    [folder,base]=OutputPath(Config);report.OutputBase=base;
    if options.pool
        exportKymo=isfield(options,'exportKymographs')&&options.exportKymographs;
        if exportKymo&&(numel(Stack)<channel||isempty(Stack{channel})||size(Stack{channel},3)<context.FrameRange(2))
            error('FIESTA:MotilityStack','Load the current-channel stack covering the requested kymograph frames.');
        end
        if numel(varargin)>1
            folder=varargin{2};
        else
            folder=uigetdir(folder,'Choose a folder for pooled motility results');if isequal(folder,0),return;end
        end
        % Each pool run has its own folder, preserving earlier exports.
        folder=fullfile(folder,[base '_motility_' char(datetime('now','Format','yyyyMMdd_HHmmssSSS'))]);mkdir(folder);
        if exportKymo
            settings=struct('Method',1,'Drift',[],'CorrectDrift',false,'StackCorrected',false,'Limits',[]);
            tools=hMainGui.RightPanel.pTools;
            if isfield(tools,'mKymoMethod'),settings.Method=get(tools.mKymoMethod,'Value');end
            if isfield(tools,'cKymoDrift'),settings.CorrectDrift=logical(get(tools.cKymoDrift,'Value'));end
            drift=getappdata(hMainGui.fig,'Drift');if numel(drift)>=channel,settings.Drift=drift{channel};end
            if isfield(hMainGui,'Menu')&&isfield(hMainGui.Menu,'mCorrectStack')
                settings.StackCorrected=strcmp(get(hMainGui.Menu.mCorrectStack,'Checked'),'on');
            end
            if isfield(hMainGui.Values,'ScaleMin'),settings.Limits=[hMainGui.Values.ScaleMin(channel) hMainGui.Values.ScaleMax(channel)];end
            report.Kymographs=fYildizMotility('savekymographs',report,hMainGui.Region,Molecule,Stack{channel},settings,folder,base);
        end
        fYildizMotility('export',report,fullfile(folder,[base '_motility.xlsx']));
        fYildizMotility('savefigures',report,folder,base);
        fprintf('Yildiz motility: %d trajectories across %d regions saved to\n%s\n',height(report.Traces),height(report.Regions),folder);
        report.OutputFolder=folder;
    end
    old=findall(0,'Tag','hMotilityGui');close(old);
    gui.fig=figure('Units','pixels','Position',[100 100 900 680],'Name','Yildiz Motility Statistics', ...
        'NumberTitle','off','MenuBar','none','ToolBar','none','Tag','hMotilityGui','Color','white');
    gui.pResults=uipanel(gui.fig,'Position',[0 .08 1 .92],'BorderType','none','BackgroundColor','white');
    gui.Report=report;gui.OutputFolder=folder;gui.OutputBase=base;
    labels={'Export XLSX','Save Figure','Close'};actions={'ExportCSV','SaveFig','Cancel'};
    for k=1:3
        action=actions{k};uicontrol(gui.fig,'Style','pushbutton','Units','normalized','Position',[.36+(k-1)*.21 .015 .19 .05], ...
            'String',labels{k},'BackgroundColor',[.16 .16 .16],'ForegroundColor',[.9 .9 .9], ...
            'Callback',@(s,e)fMotilityGui(action,getappdata(0,'hMotilityGui')));
    end
    setappdata(0,'hMotilityGui',gui);setappdata(gui.fig,'Report',report);
    setappdata(gui.fig,'CompResults',table2array(report.Traces(:,1:5)));
    setappdata(gui.fig,'ResultTitles',report.Traces.Properties.VariableNames(1:5));Draw(gui);
catch err
    fMsgDlg(err.message,'error');rethrow(err);
end

function Draw(gui)
fYildizMotility('plot',gui.Report,gui.pResults);

function ExportCSV(gui)
file=fullfile(gui.OutputFolder,[gui.OutputBase '_motility.xlsx']);
fYildizMotility('export',gui.Report,file,true);fprintf('Saved %s\n',file);

function SaveFig(gui)
fYildizMotility('savefigures',gui.Report,gui.OutputFolder,gui.OutputBase);

function [folder,base]=OutputPath(config)
folder=config.Directory;if iscell(folder),folder=folder{1};end
name=config.StackName;if iscell(name),name=name{1};end
[path,base]=fileparts(name);if ~isempty(path),folder=fullfile(folder,path);end
if ~isfolder(folder),folder=pwd;end
