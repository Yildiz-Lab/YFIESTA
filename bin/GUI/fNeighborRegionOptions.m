function opts = fNeighborRegionOptions(prevopts)
% give user GUI option to set region radii in increasing order

% c = 0.2353*ones(1,3); % for some reason usual gca doesn't work here
c = 1.2*get(gcf,'Color');
ctext = 0.9*ones(1,3);
cbutton = [0.16,0.16,0.16];
ctextbox = [0.12,0.12,0.12];
% c = get(hMainGui.fig,'Color');

hNeighborRegion.fig = figure('Units','normalized','WindowStyle','normal','DockControls','off','IntegerHandle','off','MenuBar','none','Name','Stepping Analysis Options',...
                      'NumberTitle','off','Position',[0.65 0.15 0.35 0.7],'HandleVisibility','callback','Tag','hPauseAnalysisOpts',...
                      'Visible','off','Resize','off','Color',c);
                  
fPlaceFig(hNeighborRegion.fig ,'special');

c = get(hNeighborRegion.fig ,'Color');

% make an array of grid points, say for 5 regions (though it could go on
% forever and be easily extended)

%load previous options that were saved by user
if nargin > 0
    % load function
    fprintf("Loading Previous Saved Options... \n")
    prevopts.UseNeighborRegions = 0;
else %defaults
    prevopts = struct();
    prevopts.UseNeighborRegions = 0;
    prevopts.XB = [0,100,200]; prevopts.XA = prevopts.XB;
    prevopts.YB = [0,50,100]; prevopts.YA = prevopts.YB;
    prevopts.Merge = 0; prevopts.k1 = 0; prevopts.k2=1; prevopts.FwdDwells = 0;
end
prevopts.XB(prevopts.XB == 0) = []; prevopts.XA(prevopts.XA == 0) = [];
prevopts.YB(prevopts.YB == 0) = []; prevopts.YA(prevopts.YA == 0) = [];



hNeighborRegion.cUseNeighborRegions = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.075 0.90 0.5 0.08],'Enable','on','FontSize',12,'Value',prevopts.UseNeighborRegions,...
                                          'String','Use Neighbor Regions','Style','checkbox','Tag','Neighbors','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);  

                                      
hNeighborRegion.tRegion(1) = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.2 0.785 0.12 0.08],'Enable','on','FontSize',12,...
                                 'String','Region 1','Style','text','Tag','tXB','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);                 

hNeighborRegion.tRegion(2) = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.35 0.785 0.12 0.08],'Enable','on','FontSize',12,...
                                 'String','Region 2','Style','text','Tag','tXB','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);  
               
hNeighborRegion.tRegion(3) = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.5 0.785 0.12 0.08],'Enable','on','FontSize',12,...
                                 'String','Region 3','Style','text','Tag','tXB','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);                 

hNeighborRegion.tRegion(4) = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.65 0.785 0.12 0.08],'Enable','on','FontSize',12,...
                                 'String','Region 4','Style','text','Tag','tXB','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);  

hNeighborRegion.tRegion(5) = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.8 0.785 0.12 0.08],'Enable','on','FontSize',12,...
                                 'String','Region 5','Style','text','Tag','tXB','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);  
                             

                             Z = prevopts.XB;
hNeighborRegion.tXB = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.06 0.71 0.12 0.08],'Enable','on','FontSize',12,...
                                 'String','X Before','Style','text','Tag','tXB','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);                 
                            
hNeighborRegion.eXB(1) = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.2 0.725 0.12 0.08],'Enable','on','FontSize',12,...
                                'String',UnpackageRegions(Z, 1),'Style','edit','Tag','eXB','HorizontalAlignment','center','BackgroundColor',c,'ForegroundColor',ctext);   
                            
hNeighborRegion.eXB(2) = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.35 0.725 0.12 0.08],'Enable','on','FontSize',12,...
                                'String',UnpackageRegions(Z, 2),'Style','edit','Tag','eXB','HorizontalAlignment','center','BackgroundColor',c,'ForegroundColor',ctext);   
                                                
hNeighborRegion.eXB(3) = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.5 0.725 0.12 0.08],'Enable','on','FontSize',12,...
                                'String',UnpackageRegions(Z, 3),'Style','edit','Tag','eXB','HorizontalAlignment','center','BackgroundColor',c,'ForegroundColor',ctext); 

hNeighborRegion.eXB(4) = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.65 0.725 0.12 0.08],'Enable','on','FontSize',12,...
                                'String',UnpackageRegions(Z, 4),'Style','edit','Tag','eXB','HorizontalAlignment','center','BackgroundColor',c,'ForegroundColor',ctext);   

hNeighborRegion.eXB(5) = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.8 0.725 0.12 0.08],'Enable','on','FontSize',12,...
                                'String',UnpackageRegions(Z, 5),'Style','edit','Tag','eXB','HorizontalAlignment','center','BackgroundColor',c,'ForegroundColor',ctext);   

                            Z = prevopts.XA;
hNeighborRegion.tXA = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.06 0.585 0.12 0.08],'Enable','on','FontSize',12,...
                                 'String','X After','Style','text','Tag','tXA','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);                 
                            
hNeighborRegion.eXA(1) = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.2 0.6 0.12 0.08],'Enable','on','FontSize',12,...
                                'String',UnpackageRegions(Z, 1),'Style','edit','Tag','eXA','HorizontalAlignment','center','BackgroundColor',c,'ForegroundColor',ctext);   
                            
hNeighborRegion.eXA(2) = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.35 0.6 0.12 0.08],'Enable','on','FontSize',12,...
                                'String',UnpackageRegions(Z, 2),'Style','edit','Tag','eXA','HorizontalAlignment','center','BackgroundColor',c,'ForegroundColor',ctext);   
                                                
hNeighborRegion.eXA(3) = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.5 0.6 0.12 0.08],'Enable','on','FontSize',12,...
                                'String',UnpackageRegions(Z, 3),'Style','edit','Tag','eXA','HorizontalAlignment','center','BackgroundColor',c,'ForegroundColor',ctext); 

hNeighborRegion.eXA(4) = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.65 0.6 0.12 0.08],'Enable','on','FontSize',12,...
                                'String',UnpackageRegions(Z, 4),'Style','edit','Tag','eXA','HorizontalAlignment','center','BackgroundColor',c,'ForegroundColor',ctext);   

hNeighborRegion.eXA(5) = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.8 0.6 0.12 0.08],'Enable','on','FontSize',12,...
                                'String',UnpackageRegions(Z, 5),'Style','edit','Tag','eXA','HorizontalAlignment','center','BackgroundColor',c,'ForegroundColor',ctext);   

                            Z = prevopts.YB;
hNeighborRegion.tYB = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.06 0.465 0.12 0.08],'Enable','on','FontSize',12,...
                                 'String','Y Left','Style','text','Tag','tYB','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);                 
                            
hNeighborRegion.eYB(1) = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.2 0.475 0.12 0.08],'Enable','on','FontSize',12,...
                                'String',UnpackageRegions(Z, 1),'Style','edit','Tag','eYB','HorizontalAlignment','center','BackgroundColor',c,'ForegroundColor',ctext);   
                            
hNeighborRegion.eYB(2) = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.35 0.475 0.12 0.08],'Enable','on','FontSize',12,...
                                'String',UnpackageRegions(Z, 2),'Style','edit','Tag','eYB','HorizontalAlignment','center','BackgroundColor',c,'ForegroundColor',ctext);   
                                                
hNeighborRegion.eYB(3) = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.5 0.475 0.12 0.08],'Enable','on','FontSize',12,...
                                'String',UnpackageRegions(Z, 3),'Style','edit','Tag','eYB','HorizontalAlignment','center','BackgroundColor',c,'ForegroundColor',ctext); 

hNeighborRegion.eYB(4) = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.65 0.475 0.12 0.08],'Enable','on','FontSize',12,...
                                'String',UnpackageRegions(Z, 4),'Style','edit','Tag','eYB','HorizontalAlignment','center','BackgroundColor',c,'ForegroundColor',ctext);   

hNeighborRegion.eYB(5) = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.8 0.475 0.12 0.08],'Enable','on','FontSize',12,...
                                'String',UnpackageRegions(Z, 5),'Style','edit','Tag','eYB','HorizontalAlignment','center','BackgroundColor',c,'ForegroundColor',ctext);  
                            

                            Z = prevopts.YA;
hNeighborRegion.tYA = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.06 0.34 0.12 0.08],'Enable','on','FontSize',12,...
                                 'String','Y Right','Style','text','Tag','tYA','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);                 
                            
hNeighborRegion.eYA(1) = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.2 0.35 0.12 0.08],'Enable','on','FontSize',12,...
                                'String',UnpackageRegions(Z, 1),'Style','edit','Tag','eYA','HorizontalAlignment','center','BackgroundColor',c,'ForegroundColor',ctext);   
                            
hNeighborRegion.eYA(2) = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.35 0.35 0.12 0.08],'Enable','on','FontSize',12,...
                                'String',UnpackageRegions(Z, 2),'Style','edit','Tag','eYA','HorizontalAlignment','center','BackgroundColor',c,'ForegroundColor',ctext);   
                                                
hNeighborRegion.eYA(3) = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.5 0.35 0.12 0.08],'Enable','on','FontSize',12,...
                                'String',UnpackageRegions(Z, 3),'Style','edit','Tag','eYA','HorizontalAlignment','center','BackgroundColor',c,'ForegroundColor',ctext); 

hNeighborRegion.eYA(4) = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.65 0.35 0.12 0.08],'Enable','on','FontSize',12,...
                                'String',UnpackageRegions(Z, 4),'Style','edit','Tag','eYA','HorizontalAlignment','center','BackgroundColor',c,'ForegroundColor',ctext);   

hNeighborRegion.eYA(5) = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.8 0.35 0.12 0.08],'Enable','on','FontSize',12,...
                                'String',UnpackageRegions(Z, 5),'Style','edit','Tag','eYA','HorizontalAlignment','center','BackgroundColor',c,'ForegroundColor',ctext);  

                          
                            % JS Edit 2024/05/14 giving options for fits

% hNeighborRegion.tMerge = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.05 0.3 0.3 0.05],'Enable','on','FontSize',12,...
%                                  'String','Merge on-off axis steps','Style','text','Tag','tYA','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);                
                             
hNeighborRegion.cMerge = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.1 0.3 0.3 0.05],'Enable','on','FontSize',12,'Value',prevopts.Merge,...
                                          'String','Merge on-off axis steps','Style','checkbox','Tag','Neighbors','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);  


hNeighborRegion.tLifetime = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.05 0.25 0.3 0.05],'Enable','on','FontSize',12,...
                                 'String','Rate fitting options','Style','text','Tag','tYA','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);                
                             
hNeighborRegion.cPoissonk1 = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.1 0.2 0.3 0.05],'Enable','on','FontSize',12,'Value',prevopts.k1,...
                                          'String','Exponential Decay','Style','checkbox','Tag','Neighbors','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);  

hNeighborRegion.cPoissonk2 = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.4 0.2 0.3 0.05],'Enable','on','FontSize',12,'Value',prevopts.k2,...
                                          'String','Poisson k=2','Style','checkbox','Tag','Neighbors','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext); 

hNeighborRegion.cFwdDwells = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.63 0.2 0.37 0.05],'Enable','on','FontSize',12,'Value',prevopts.FwdDwells,...
                                          'String','Use forward dwells only','Style','checkbox','Tag','Neighbors','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext); 

                            % Buttons
hNeighborRegion.bOkay = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.55 0.05 0.42 0.1],'Enable','on','FontSize',14,...
                                'String','Okay','Style','pushbutton','Tag','bOkay','HorizontalAlignment','center','BackgroundColor',cbutton,'ForegroundColor',ctext,...
                                'Callback',@Okay);  
                            
hNeighborRegion.bSetParams = uicontrol('Parent',hNeighborRegion.fig,'Units','normalized','Position',[0.1 0.05 0.42 0.1],'Enable','on','FontSize',14,...
                                'String','Set Regions & Proceed','Style','pushbutton','Tag','bFindThresh','HorizontalAlignment','center','BackgroundColor',cbutton,'ForegroundColor',ctext,...
                                'Callback',@SetRegions);  


setappdata(0,'hNeighborRegionOpts',hNeighborRegion);
setappdata(hNeighborRegion.fig,'opts',[]);
uiwait(hNeighborRegion.fig)
try
    opts = getappdata(hNeighborRegion.fig,'opts');
    close(hNeighborRegion.fig);
catch 
    opts = [];
end


function SetRegions(~,~)
global Config;
hNeighborRegion = getappdata(0,'hNeighborRegionOpts');
opts = PackageRegions();
opts.Merge = hNeighborRegion.cMerge.Value;
opts.Poissonk1 = hNeighborRegion.cPoissonk1.Value;
opts.Poissonk2 = hNeighborRegion.cPoissonk2.Value;
opts.FwdDwells = hNeighborRegion.cFwdDwells.Value;
Config.NeighborRegionsOpts = opts;
opts.UseNeighborRegions = hNeighborRegion.cUseNeighborRegions.Value;
setappdata(hNeighborRegion.fig,'opts',opts);
uiresume(gcbf);


function Okay(~,~)
hNeighborRegion = getappdata(0,'hNeighborRegionOpts');
opts = PackageRegions();
opts.UseNeighborRegions = hNeighborRegion.cUseNeighborRegions.Value;
opts.Merge = hNeighborRegion.cMerge.Value;
opts.Poissonk1 = hNeighborRegion.cPoissonk1.Value;
opts.Poissonk2 = hNeighborRegion.cPoissonk2.Value;
opts.FwdDwells = hNeighborRegion.cFwdDwells.Value;
setappdata(hNeighborRegion.fig,'opts',opts);
uiresume(gcbf);


function regions = PackageRegions()

hNeighborRegion = getappdata(0,'hNeighborRegionOpts');

stopidx = 0;
for i = 1:length(hNeighborRegion.eXB)
    if (~isempty(hNeighborRegion.eXB(i).String)) || (~isempty(hNeighborRegion.eXA(i).String)) || (~isempty(hNeighborRegion.eYB(i).String)) || (~isempty(hNeighborRegion.eYA(i).String))
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
    if isempty(str2num(hNeighborRegion.eXB(j).String))
        XB(j+1) = 0;
    else
        XB(j+1) = str2num(hNeighborRegion.eXB(j).String);
    end
    if XB(j+1) < XB(j)
        XB(j+1) = XB(j);
    end
    
    if isempty(str2num(hNeighborRegion.eXA(j).String))
        XA(j+1) = 0;
    else
        XA(j+1) = str2num(hNeighborRegion.eXA(j).String);
    end
    if XA(j+1) < XA(j)
        XA(j+1) = XA(j);
    end
    
    if isempty(str2num(hNeighborRegion.eYB(j).String))
        YB(j+1) = 0;
    else
        YB(j+1) = str2num(hNeighborRegion.eYB(j).String);
    end
    if YB(j+1) < YB(j)
        YB(j+1) = YB(j);
    end
    
    if isempty(str2num(hNeighborRegion.eYA(j).String))
        YA(j+1) = 0;
    else
        YA(j+1) = str2num(hNeighborRegion.eYA(j).String);
    end
    if YA(j+1) < YA(j)
        YA(j+1) = YA(j);
    end

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






