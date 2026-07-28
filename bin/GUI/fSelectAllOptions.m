function choice = fSelectAllOptions
% give user GUI option to set Neighbor step finding options if selected
% Adapted from fNeighborRegionOptions

% c = 0.2353*ones(1,3); % for some reason usual gca doesn't work here
c = 1.2*get(gcf,'Color');
ctext = 0.9*ones(1,3);
cbutton = [0.16,0.16,0.16];
ctextbox = [0.12,0.12,0.12];
% c = get(hMainGui.fig,'Color');

hSelectOpts.fig = figure('Units','normalized','WindowStyle','normal','DockControls','off','IntegerHandle','off','MenuBar','none','Name','Select Options',...
                      'NumberTitle','off','Position',[0.65 0.15 0.35 0.7],'HandleVisibility','callback','Tag','hPauseAnalysisOpts',...
                      'Visible','off','Resize','off','Color',c);
                  
fPlaceFig(hSelectOpts.fig ,'special');

c = get(hSelectOpts.fig ,'Color');

                            % Buttons to proceed
hSelectOpts.bAll = uicontrol('Parent',hSelectOpts.fig,'Units','normalized','Position',[0.2 0.65 0.62 0.1],'Enable','on','FontSize',14,...
                                'String','Select all','Style','pushbutton','Tag','bPointSelect','HorizontalAlignment','center','BackgroundColor',cbutton,'ForegroundColor',ctext,...
                                'Callback',@fPointSelectAll);  
                            
hSelectOpts.bBefore = uicontrol('Parent',hSelectOpts.fig,'Units','normalized','Position',[0.2 0.45 0.62 0.1],'Enable','on','FontSize',14,...
                                'String','Select all before','Style','pushbutton','Tag','bPointSelect','HorizontalAlignment','center','BackgroundColor',cbutton,'ForegroundColor',ctext,...
                                'Callback',@fPointSelectBefore); 

hSelectOpts.bAfter = uicontrol('Parent',hSelectOpts.fig,'Units','normalized','Position',[0.2 0.25 0.62 0.1],'Enable','on','FontSize',14,...
                                'String','Select all after','Style','pushbutton','Tag','bPointSelect','HorizontalAlignment','center','BackgroundColor',cbutton,'ForegroundColor',ctext,...
                                'Callback',@fPointSelectAfter); 


setappdata(0,'hSelectOpts',hSelectOpts);
setappdata(hSelectOpts.fig,'choice',0);
uiwait(hSelectOpts.fig)
try
    choice = getappdata(hSelectOpts.fig,'choice');
    close(hSelectOpts.fig);
catch 
    choice = 0;
end


function fPointSelectAll(~,~)
hSelectOpts = getappdata(0,'hSelectOpts');
setappdata(hSelectOpts.fig,'choice',0);
uiresume(gcbf);

function fPointSelectBefore(~,~)
hSelectOpts = getappdata(0,'hSelectOpts');
setappdata(hSelectOpts.fig,'choice',1);
uiresume(gcbf);

function fPointSelectAfter(~,~)
hSelectOpts = getappdata(0,'hSelectOpts');
setappdata(hSelectOpts.fig,'choice',2);
uiresume(gcbf);
