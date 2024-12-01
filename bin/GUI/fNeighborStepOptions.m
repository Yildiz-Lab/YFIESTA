function opts = fNeighborStepOptions(prevopts)
% give user GUI option to set Neighbor step finding options if selected
% Adapted from fNeighborRegionOptions

% c = 0.2353*ones(1,3); % for some reason usual gca doesn't work here
c = 1.2*get(gcf,'Color');
ctext = 0.9*ones(1,3);
cbutton = [0.16,0.16,0.16];
ctextbox = [0.12,0.12,0.12];
% c = get(hMainGui.fig,'Color');

hNeighborStep.fig = figure('Units','normalized','WindowStyle','normal','DockControls','off','IntegerHandle','off','MenuBar','none','Name','Stepping Analysis Options',...
                      'NumberTitle','off','Position',[0.65 0.15 0.35 0.7],'HandleVisibility','callback','Tag','hPauseAnalysisOpts',...
                      'Visible','off','Resize','off','Color',c);
                  
fPlaceFig(hNeighborStep.fig ,'special');

c = get(hNeighborStep.fig ,'Color');

% make an array of grid points, say for 5 regions (though it could go on
% forever and be easily extended)

%load previous options that were saved by user
if nargin > 0
    % load function
    fprintf("Loading Previous Saved Options... \n")
else %defaults
    prevopts = struct();
    prevopts.XB = [200]; prevopts.XA = prevopts.XB;
    % prevopts.YB = [5]; prevopts.YA = prevopts.YB;
    prevopts.ExcludeTime = 0;
    prevopts.eExcludeTime = num2str('0');
    prevopts.ExistThresh = num2str('0.4');
    prevopts.InterpRes = num2str('0.25');
    prevopts.ShowPathExt = 0;
end
prevopts.XB(prevopts.XB == 0) = []; prevopts.XA(prevopts.XA == 0) = [];
% prevopts.YB(prevopts.YB == 0) = []; prevopts.YA(prevopts.YA == 0) = [];


                             

                             Z = prevopts.XB;
hNeighborStep.tXB = uicontrol('Parent',hNeighborStep.fig,'Units','normalized','Position',[0.06 0.85 0.12 0.08],'Enable','on','FontSize',12,...
                                 'String','X Before','Style','text','Tag','tXB','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);                 
                            
hNeighborStep.eXB(1) = uicontrol('Parent',hNeighborStep.fig,'Units','normalized','Position',[0.2 0.865 0.12 0.08],'Enable','on','FontSize',12,...
                                'String',UnpackageRegions(Z, 1),'Style','edit','Tag','eXB','HorizontalAlignment','center','BackgroundColor',c,'ForegroundColor',ctext);   
                            
% hNeighborStep.eXB(2) = uicontrol('Parent',hNeighborStep.fig,'Units','normalized','Position',[0.35 0.765 0.12 0.08],'Enable','on','FontSize',12,...
%                                 'String',UnpackageRegions(Z, 2),'Style','edit','Tag','eXB','HorizontalAlignment','center','BackgroundColor',c,'ForegroundColor',ctext);   

                            Z = prevopts.XA;
hNeighborStep.tXA = uicontrol('Parent',hNeighborStep.fig,'Units','normalized','Position',[0.4 0.85 0.12 0.08],'Enable','on','FontSize',12,...
                                 'String','X After','Style','text','Tag','tXA','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);                 
                            
hNeighborStep.eXA(1) = uicontrol('Parent',hNeighborStep.fig,'Units','normalized','Position',[0.54 0.865 0.12 0.08],'Enable','on','FontSize',12,...
                                'String',UnpackageRegions(Z, 1),'Style','edit','Tag','eXA','HorizontalAlignment','center','BackgroundColor',c,'ForegroundColor',ctext);   
                            
% hNeighborStep.eXA(2) = uicontrol('Parent',hNeighborStep.fig,'Units','normalized','Position',[0.35 0.64 0.12 0.08],'Enable','on','FontSize',12,...
%                                 'String',UnpackageRegions(Z, 2),'Style','edit','Tag','eXA','HorizontalAlignment','center','BackgroundColor',c,'ForegroundColor',ctext);   

%                             Z = prevopts.YB;
% hNeighborStep.tYB = uicontrol('Parent',hNeighborStep.fig,'Units','normalized','Position',[0.06 0.505 0.12 0.08],'Enable','on','FontSize',12,...
%                                  'String','Y Left','Style','text','Tag','tYB','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);                 
% 
% hNeighborStep.eYB(1) = uicontrol('Parent',hNeighborStep.fig,'Units','normalized','Position',[0.2 0.515 0.12 0.08],'Enable','on','FontSize',12,...
%                                 'String',UnpackageRegions(Z, 1),'Style','edit','Tag','eYB','HorizontalAlignment','center','BackgroundColor',c,'ForegroundColor',ctext);   
                            
% hNeighborStep.eYB(2) = uicontrol('Parent',hNeighborStep.fig,'Units','normalized','Position',[0.35 0.515 0.12 0.08],'Enable','on','FontSize',12,...
%                                 'String',UnpackageRegions(Z, 2),'Style','edit','Tag','eYB','HorizontalAlignment','center','BackgroundColor',c,'ForegroundColor',ctext);   
                            

%                             Z = prevopts.YA;
% hNeighborStep.tYA = uicontrol('Parent',hNeighborStep.fig,'Units','normalized','Position',[0.06 0.38 0.12 0.08],'Enable','on','FontSize',12,...
%                                  'String','Y Right','Style','text','Tag','tYA','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);                 
% 
% hNeighborStep.eYA(1) = uicontrol('Parent',hNeighborStep.fig,'Units','normalized','Position',[0.2 0.39 0.12 0.08],'Enable','on','FontSize',12,...
%                                 'String',UnpackageRegions(Z, 1),'Style','edit','Tag','eYA','HorizontalAlignment','center','BackgroundColor',c,'ForegroundColor',ctext);   
                            
% hNeighborStep.eYA(2) = uicontrol('Parent',hNeighborStep.fig,'Units','normalized','Position',[0.35 0.39 0.12 0.08],'Enable','on','FontSize',12,...
%                                 'String',UnpackageRegions(Z, 2),'Style','edit','Tag','eYA','HorizontalAlignment','center','BackgroundColor',c,'ForegroundColor',ctext);   

                          
                            % JS Edit 2024/11/30 giving options for step display of neighbors

hNeighborStep.tInterpRes = uicontrol('Parent',hNeighborStep.fig,'Units','normalized','Position',[0.1 0.765 0.45 0.05],'Enable','on','FontSize',12,...
                                          'String','Interpolation Spacing (nm)','Style','text','Tag','Neighbors','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);  
hNeighborStep.eInterpRes = uicontrol('Parent',hNeighborStep.fig,'Units','normalized','Position',[0.5 0.75 0.12 0.08],'Enable','on','FontSize',12,'String',prevopts.InterpRes,...
                                'Style','edit','Tag','eYA','HorizontalAlignment','center','BackgroundColor',c,'ForegroundColor',ctext);

hNeighborStep.tExcludeTime = uicontrol('Parent',hNeighborStep.fig,'Units','normalized','Position',[0.1 0.665 0.35 0.05],'Enable','on','FontSize',12,...
                                          'String','Neighbor Exist Thresh (s)','Style','text','Tag','Neighbors','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);  
hNeighborStep.eExistThresh = uicontrol('Parent',hNeighborStep.fig,'Units','normalized','Position',[0.5 0.65 0.12 0.08],'Enable','on','FontSize',12,'String',prevopts.ExistThresh,...
                                'Style','edit','Tag','eYA','HorizontalAlignment','center','BackgroundColor',c,'ForegroundColor',ctext);

hNeighborStep.cExcludeTime = uicontrol('Parent',hNeighborStep.fig,'Units','normalized','Position',[0.1 0.565 0.35 0.05],'Enable','on','FontSize',12,'Value',prevopts.ExcludeTime,...
                                          'String','Exclude by time (s)','Style','checkbox','Tag','Neighbors','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);
hNeighborStep.eExcludeTime = uicontrol('Parent',hNeighborStep.fig,'Units','normalized','Position',[0.5 0.55 0.12 0.08],'Enable','on','FontSize',12,'String',prevopts.eExcludeTime,...
                                'Style','edit','Tag','eYA','HorizontalAlignment','center','BackgroundColor',c,'ForegroundColor',ctext);


hNeighborStep.cShowPathExt = uicontrol('Parent',hNeighborStep.fig,'Units','normalized','Position',[0.1 0.465 0.35 0.05],'Enable','on','FontSize',12,'Value',prevopts.ShowPathExt,...
                                          'String','Show Path Extension','Style','checkbox','Tag','Neighbors','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);



                            % Buttons to proceed
hNeighborStep.bOkay = uicontrol('Parent',hNeighborStep.fig,'Units','normalized','Position',[0.55 0.05 0.42 0.1],'Enable','on','FontSize',14,...
                                'String','Okay','Style','pushbutton','Tag','bOkay','HorizontalAlignment','center','BackgroundColor',cbutton,'ForegroundColor',ctext,...
                                'Callback',@Okay);  
                            
hNeighborStep.bSetParams = uicontrol('Parent',hNeighborStep.fig,'Units','normalized','Position',[0.1 0.05 0.42 0.1],'Enable','on','FontSize',14,...
                                'String','Set Opts & Proceed','Style','pushbutton','Tag','bSetParams','HorizontalAlignment','center','BackgroundColor',cbutton,'ForegroundColor',ctext,...
                                'Callback',@SetRegions);  


setappdata(0,'hNeighborStepOpts',hNeighborStep);
setappdata(hNeighborStep.fig,'opts',[]);
uiwait(hNeighborStep.fig)
try
    opts = getappdata(hNeighborStep.fig,'opts');
    close(hNeighborStep.fig);
catch 
    opts = [];
end


function SetRegions(~,~)
global Config;
hNeighborStep = getappdata(0,'hNeighborStepOpts');
opts = PackageRegions();
opts.ExcludeTime = hNeighborStep.cExcludeTime.Value;
opts.eExcludeTime = str2num(hNeighborStep.eExcludeTime.String);
opts.ExistThresh = str2num(hNeighborStep.eExistThresh.String);
opts.InterpRes = str2num(hNeighborStep.eInterpRes.String);
opts.ShowPathExt = hNeighborStep.cShowPathExt.Value;
Config.NeighborStepsOpts = opts;
setappdata(hNeighborStep.fig,'opts',opts);
uiresume(gcbf);


function Okay(~,~)
hNeighborStep = getappdata(0,'hNeighborStepOpts');
opts = PackageRegions();
opts.ExcludeTime = hNeighborStep.cExcludeTime.Value;
opts.eExcludeTime = str2num(hNeighborStep.eExcludeTime.String);
opts.ExistThresh = str2num(hNeighborStep.eExistThresh.String);
opts.InterpRes = str2num(hNeighborStep.eInterpRes.String);
opts.ShowPathExt = hNeighborStep.cShowPathExt.Value;
setappdata(hNeighborStep.fig,'opts',opts);
uiresume(gcbf);


function regions = PackageRegions()

hNeighborStep = getappdata(0,'hNeighborStepOpts');

stopidx = 0;
for i = 1:length(hNeighborStep.eXB)
    if (~isempty(hNeighborStep.eXB(i).String)) || (~isempty(hNeighborStep.eXA(i).String)) || (~isempty(hNeighborStep.eYB(i).String)) || (~isempty(hNeighborStep.eYA(i).String))
        stopidx = i; %we want the highest nonempty one
    end
end

% Now we know how many regions we have to go and package into arrays.
% The rules are:
%   - If after is empty when your corresponding before is filled in, it
%   becomes that number
%   - If your number is lower than or empty compared to the number previous,
%   you make it equivalent to the previous
XB = zeros(1,stopidx+1); XA = zeros(1,stopidx+1); 
YB = zeros(1,stopidx+1); YA = zeros(1,stopidx+1); 
for j = 1:stopidx % the first should always be zero
    % do this for each dimension
    if isempty(str2num(hNeighborStep.eXB(j).String))
        XB(j+1) = 0;
    else
        XB(j+1) = str2num(hNeighborStep.eXB(j).String);
    end
    if XB(j+1) < XB(j)
        XB(j+1) = XB(j);
    end

    if isempty(str2num(hNeighborStep.eXA(j).String))
        XA(j+1) = 0;
    else
        XA(j+1) = str2num(hNeighborStep.eXA(j).String);
    end
    if XA(j+1) < XA(j)
        XA(j+1) = XA(j);
    end

    % if isempty(str2num(hNeighborStep.eYB(j).String))
    %     YB(j+1) = 0;
    % else
    %     YB(j+1) = str2num(hNeighborStep.eYB(j).String);
    % end
    % if YB(j+1) < YB(j)
    %     YB(j+1) = YB(j);
    % end
    % 
    % if isempty(str2num(hNeighborStep.eYA(j).String))
    %     YA(j+1) = 0;
    % else
    %     YA(j+1) = str2num(hNeighborStep.eYA(j).String);
    % end
    % if YA(j+1) < YA(j)
    %     YA(j+1) = YA(j);
    % end

end

% now XB, XA, YB, YA have their normal meanings which means we just have to
% pass them into fNeighborlyRegions eventually. Hence, we'll just save in
% parameters
regions.XB = XB; regions.YB = YB; regions.XA = XA; regions.YA = YA;


function fillstr = UnpackageRegions(Z, i)
% Z : XB, XA, YB, YA
% i : region i you want to fill
stopidx = length(Z);
if i > stopidx
    fillstr = '';
else
    fillstr = num2str(Z(i));
end






