function opts = fPauseAnalysisOpts(prevopts)
% JS 2023/01/27
% options for pause analysis gui

c = 1.2*get(gcf,'Color');
ctext = 0.9*ones(1,3);
cbutton = [0.16,0.16,0.16];
ctextbox = [0.12,0.12,0.12];
%c = get(hMainGui.fig,'Color');

hPauseAnalysisOpts.fig = figure('Units','normalized','WindowStyle','normal','DockControls','off','IntegerHandle','off','MenuBar','none','Name','Pause Analysis Options',...
                      'NumberTitle','off','Position',[0.65 0.15 0.35 0.7],'HandleVisibility','callback','Tag','hPauseAnalysisOpts',...
                      'Visible','off','Resize','off','Color',c);
                  
fPlaceFig(hPauseAnalysisOpts.fig ,'small');

c = get(hPauseAnalysisOpts.fig ,'Color');

%load previous options that were saved by user
if nargin > 0
    % load function
    fprintf("Loading Previous Saved Options... \n")
    prevopts.UseNeighborRegions = 0;
else %defaults
    prevopts = struct();
    prevopts.UseNeighborRegions = 0;
    prevopts.PauseThreshold = 47;
end


hPauseAnalysisOpts.rSetPauseThreshtext = uicontrol('Parent',hPauseAnalysisOpts.fig,'Units','normalized','Position',[0.1 0.715 0.45 0.2],'Enable','on','FontSize',12,...
                                        'String','Pause Event Threshold','Style','text','Tag','SetPauseThresh','BackgroundColor',c,'ForegroundColor',ctext);

hPauseAnalysisOpts.rSetPauseThresh = uicontrol('Parent',hPauseAnalysisOpts.fig,'Units','normalized','Position',[0.6 0.725 0.2 0.2],'Enable','on','FontSize',12,...
                                        'String',num2str(prevopts.PauseThreshold),'Style','edit','Tag','SetPauseThresh','BackgroundColor',cbutton,'ForegroundColor',ctext);

hPauseAnalysisOpts.cUseNeighborRegions = uicontrol('Parent',hPauseAnalysisOpts.fig,'Units','normalized','Position',[0.1 0.45 0.9 0.25],'Enable','on','FontSize',12,'Value',prevopts.UseNeighborRegions,...
                                          'String','Use Neighbor Regions','Style','checkbox','Tag','Neighbors','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);  
    
hPauseAnalysisOpts.bOkay = uicontrol('Parent',hPauseAnalysisOpts.fig,'Units','normalized','Position',[0.55 0.1 0.35 0.2],'Enable','on','FontSize',14,...
                                'String','Okay','Style','pushbutton','Tag','bOkay','HorizontalAlignment','center','BackgroundColor',cbutton,'ForegroundColor',ctext,...
                                'Callback',@Okay);  
                            
hPauseAnalysisOpts.bFindThresh = uicontrol('Parent',hPauseAnalysisOpts.fig,'Units','normalized','Position',[0.1 0.1 0.35 0.2],'Enable','on','FontSize',14,...
                                'String','Set Reference Threshold','Style','pushbutton','Tag','bFindThresh','HorizontalAlignment','center','BackgroundColor',cbutton,'ForegroundColor',ctext,...
                                'Callback',@FindThresh);  
                            
                           
setappdata(0,'hPauseAnalysisOpts',hPauseAnalysisOpts);
setappdata(hPauseAnalysisOpts.fig,'opts',[]);
uiwait(hPauseAnalysisOpts.fig)
try
    opts = getappdata(hPauseAnalysisOpts.fig,'opts');
    close(hPauseAnalysisOpts.fig);
catch 
    opts =[];
end


function FindThresh(~,~)
hPauseAnalysisOpts = getappdata(0,'hPauseAnalysisOpts');
opts.UseNeighborRegions = hPauseAnalysisOpts.cUseNeighborRegions.Value;
opts.PauseThresh = 0;
setappdata(hPauseAnalysisOpts.fig,'opts',opts);
uiresume(gcbf);


function Okay(~,~)
hPauseAnalysisOpts = getappdata(0,'hPauseAnalysisOpts');
opts.UseNeighborRegions = hPauseAnalysisOpts.cUseNeighborRegions.Value;
opts.PauseThresh = str2double(hPauseAnalysisOpts.rSetPauseThresh.String);
setappdata(hPauseAnalysisOpts.fig,'opts',opts);
uiresume(gcbf);
