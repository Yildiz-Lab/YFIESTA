function options = fMotilityStatsDlg
% JS Function 2022/07/14 to pass in options for yildiz motility assays


hMotilityStatsDlg = dialog('Name','Options for FIESTA Yildiz Motility Statistics','Visible','off');
fPlaceFig(hMotilityStatsDlg,'speed');


h(1) = uicontrol('Parent',hMotilityStatsDlg ,'Units','normalized','Position',[0.05 0.92 0.9 0.08],'Style','text','Tag','tFit',...
          'String','Thresholds for a moving molecule','HorizontalAlignment','left','FontSize',12);
              
h(2) = uicontrol('Parent',hMotilityStatsDlg ,'Units','normalized','Position',[0.05 0.85 0.7 0.08],'Style','text','Tag','tFit',...
          'String','Distance Threshold (pixels)','HorizontalAlignment','left','FontSize',10);
      
h(3) = uicontrol('Parent',hMotilityStatsDlg ,'Units','normalized','Position',[0.75 0.85 0.2 0.08],'Style','edit','Tag','edist',...
          'String','3','HorizontalAlignment','left','FontSize',10);

h(4) = uicontrol('Parent',hMotilityStatsDlg ,'Units','normalized','Position',[0.05 0.75 0.9 0.08],'Style','text','Tag','tFit',...
          'String','Threshold for molecule existence','HorizontalAlignment','left','FontSize',12);
      
h(5) = uicontrol('Parent',hMotilityStatsDlg ,'Units','normalized','Position',[0.05 0.68 0.7 0.08],'Style','text','Tag','tFit',...
          'String','Time Threshold (pixels)','HorizontalAlignment','left','FontSize',10);
      
h(6) = uicontrol('Parent',hMotilityStatsDlg ,'Units','normalized','Position',[0.75 0.68 0.2 0.08],'Style','edit','Tag','edist',...
          'String','2','HorizontalAlignment','left','FontSize',10);
      
h(7) = uicontrol('Parent',hMotilityStatsDlg ,'Units','normalized','Position',[0.05 0.58 0.9 0.08],'Style','text','Tag','tFit',...
          'String','Truncate ends (NOT FUNCTIONAL)','HorizontalAlignment','left','FontSize',12);
      
h(8) = uicontrol('Parent',hMotilityStatsDlg ,'Units','normalized','Position',[0.05 0.51 0.7 0.08],'Style','text','Tag','tFit',...
          'String','Distance Threshold (pixels)','HorizontalAlignment','left','FontSize',10);
      
h(9) = uicontrol('Parent',hMotilityStatsDlg ,'Units','normalized','Position',[0.75 0.51 0.2 0.08],'Style','edit','Tag','edist',...
          'String','5','HorizontalAlignment','left','FontSize',10);
      
h(10) = uicontrol('Parent',hMotilityStatsDlg ,'Units','normalized','Position',[0.05 0.41 0.9 0.08],'Style','checkbox','Tag','cVelMove',...
          'String','Use only movers to calculate velocity?','HorizontalAlignment','left','FontSize',10,'Visible','on','Value',1);
      
h(11) = uicontrol('Parent',hMotilityStatsDlg ,'Units','normalized','Position',[0.05 0.33 0.9 0.08],'Style','checkbox','Tag','cVelMove',...
          'String','Use only movers to calculate run length?','HorizontalAlignment','left','FontSize',10,'Visible','on','Value',1);
      
h(12) = uicontrol('Parent',hMotilityStatsDlg ,'Units','normalized','Position',[0.05 0.25 0.9 0.08],'Style','checkbox','Tag','cVelMove',...
          'String','Use only movers to calculate landing rate (Run Flux)?','HorizontalAlignment','left','FontSize',10,'Visible','on');
      
h(13) = uicontrol('Parent',hMotilityStatsDlg ,'Units','normalized','Position',[0.05 0.17 0.9 0.08],'Style','checkbox','Tag','cVelMove',...
          'String','Ignore Edge Molecules (NOT FUNCTIONAL)?','HorizontalAlignment','left','FontSize',10,'Visible','on');


uicontrol('Parent',hMotilityStatsDlg ,'Units','normalized','Position',[0.05 0.05 0.4 0.1],'Style','pushbutton','String','Ok','FontSize',12,'Callback',@doControlCallback);
uicontrol('Parent',hMotilityStatsDlg ,'Units','normalized','Position',[0.55 0.05 0.4 0.1],'Style','pushbutton','String','Cancel','FontSize',12,'Callback',@doControlCallback);

uiwait(hMotilityStatsDlg);
if ~ishandle(hMotilityStatsDlg)
    options = [];
else
    button = get(hMotilityStatsDlg,'UserData');
    if strcmp(button,'Ok') 
        options.dthresh = str2double(get(h(3), 'String'));
        options.fthresh = str2double(get(h(6), 'String'));
        options.edgeregion = str2double(get(h(9), 'String'));
        options.velocity_moverstoggle = get(h(10), 'Value');
        options.runlength_moverstoggle = get(h(11), 'Value');
        options.landingrate_moverstoggle = get(h(12), 'Value');
        options.runlength_countedges = get(h(13), 'Value');
    else
        options = [];
    end
    delete(hMotilityStatsDlg);
end

function Update(obj, ~)
mode = obj.Value;
h = obj.UserData;
set(h,'Visible','off');
if mode == 1
    set(h(1:2),'Visible','on');
elseif mode == 2
    set(h(3:4),'Visible','on');
else
    set(h(5:6),'Visible','on');
end
 
function doControlCallback(obj, ~)
set(gcbf,'UserData',get(obj,'String'));
uiresume(gcbf);