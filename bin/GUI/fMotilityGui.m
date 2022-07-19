function fMotilityGui(func,varargin)
% JS 2022/07/14 For Yildiz Motility Assays

switch func
    case 'Create'
        Create;       
    case 'Draw'
        Draw(varargin{1});  
    case 'SaveFig'
        SaveFig(varargin{1}); 
    case 'ExportCSV'
        ExportCSV(varargin{1});
    case 'Cancel'
        Cancel(varargin{1});          
end

function Create
global Molecule;
global Filament;
global Index;
global Config;

global KymoTrackMol; %I think this may allow for selecting the tracks in a kymograph but unsure yet

h=findobj('Tag','hMotilityGui');
close(h)

hMainGui=getappdata(0,'hMainGui');

MolSelect = [Molecule.Selected];
FilSelect = [Filament.Selected];
if all(MolSelect==0) && all(FilSelect==0)
    fMsgDlg('No track selected!','error');
    return;
end

% choose options for what to include
options = fMotilityStatsDlg;

if isempty(options)
    return
end
% if strcmp(options.mode,'average')
%     AverageDis = options.dis;
% else
%     AverageDis  = 0;
% end


% Draw the figure

hMotilityGui.fig = figure('Units','normalized','DockControls','off','IntegerHandle','off','MenuBar','none','Name','Motility Statistics',...
                      'NumberTitle','off','HandleVisibility','callback','Tag','hMotilityGui',...
                      'Visible','off','Resize','off','WindowStyle','modal');
                  
if ispc
    set(hMotilityGui.fig,'Color',[236 233 216]/255);
end

c=get(hMotilityGui.fig,'Color');

hMotilityGui.pResults = uipanel('Parent',hMotilityGui.fig,'Position',[0.05 0.05 0.9 0.9],'Tag','PlotPanel','BackgroundColor','white');

hMotilityGui.aPlot = axes('Parent',hMotilityGui.pResults,'Units','normalized','OuterPosition',[0 0 1 1],'Tag','aPlotXYZ','Visible','off'); 
                         
hMotilityGui.bExportCSV = uicontrol('Parent',hMotilityGui.fig,'Units','normalized','Callback','fMotilityGui(''ExportCSV'',getappdata(0,''hMotilityGui''));',...
                            'Position',[0.4 0.01 0.175 0.03],'String','Export','Tag','bExportCSV'); 
                         
hMotilityGui.bSaveFig = uicontrol('Parent',hMotilityGui.fig,'Units','normalized','Callback','fMotilityGui(''SaveFig'',getappdata(0,''hMotilityGui''));',...
                             'Position',[0.6 0.01 0.175 0.03],'String','Save Figure','Tag','bSaveFig'); 

hMotilityGui.bCancel = uicontrol('Parent',hMotilityGui.fig,'Units','normalized','Callback','fMotilityGui(''Cancel'',getappdata(0,''hMotilityGui''));',...
                             'Position',[0.8 0.01 0.175 0.03],'String','Cancel','Tag','bCancel');


setappdata(0,'hMotilityGui',hMotilityGui);
                         
Index = find(MolSelect==1);

% Do the analysis on selected molecules

% Results: Molecule Index, Velocity, Run Length, Landing Rate
Results = [];
% this is where thresholding occurs by options
for i = Index
    %frame, time(s), x(nm), y(nm), z(nm), dist to origin(nm), FWHM(nm), Amplitude(cnts), Position Error(nm), Tags
    rr = Molecule(i).Results;
    
    % filter if they are over the dthresh and fthresh
    if rr(end,6)-rr(1,6) < options.dthresh*Config.PixSize
        moverthresh = 0;
    else
        moverthresh = 1;
    end
    
    if rr(end,1)-rr(1,1) < options.fthresh
        molthresh = 0;
    else
        molthresh = 1;
    end
    
    % then add to the associated results iff they are above the passthresh
    % for moving
    Results = [Results; nan(1,4)]; Results(end,1)=i;
    ResultTitles = ["Molecule", "Velocity (nm/s)", "Distance (nm)", "Landing Rate (um^{-1} s^{-1})"];
    if (~options.velocity_moverstoggle || moverthresh) && molthresh
    Results(end,2) = (rr(end,6) - rr(1,6))/(rr(end,2) - rr(1,2));
    end
    if (~options.runlength_moverstoggle || moverthresh) && molthresh
    Results(end,3) = rr(end,6) - rr(1,6);
    end
    if (~options.landingrate_moverstoggle || moverthresh) && molthresh
    Results(end,4) = 1;
    end
    
end

% at the end, do modification for the Landing Rate
deltad = hMainGui.Scan.InterpD(end); % distance in pixels
%distance in time +1 since the number of frames is actually the difference
deltat = str2double(get(hMainGui.RightPanel.pTools.eKymoEnd,'String')) - str2double(get(hMainGui.RightPanel.pTools.eKymoStart,'String')) + 1;

% only fill in results if not nan
idxnan = find(isnan(Results(:,4)));
if length(Config.Time) > 1
    Results(:,4) = ones(size(Results(:,1),1),1) * sum(Results(:,4),'omitnan') / deltad / Config.PixSize / deltat / Config.Time(1) * 1e6;
else
    Results(:,4) = ones(size(Results(:,1),1),1) * sum(Results(:,4),'omitnan') / deltad / Config.PixSize / deltat / Config.Time * 1e6;
end
Results(idxnan,4) = NaN;

% set this data for later
setappdata(hMotilityGui.fig,'CompResults',Results);
setappdata(hMotilityGui.fig,'ResultTitles',ResultTitles);
Draw(hMotilityGui);  

fPlaceFig(hMotilityGui.fig,'big');

function Draw(hMotilityGui)

set(hMotilityGui.aPlot,'Visible','on');   
Results = getappdata(hMotilityGui.fig,'CompResults')
Titles = getappdata(hMotilityGui.fig,'ResultTitles');


%plots
for i=2:size(Results,2)
    subplot(1,size(Results,2)-1,i-1);
    cdfplot(Results(:,i));
    title(Titles(i))
    legend("Mean: " +num2str(mean(Results(:,i),'omitnan')) + newline + "Std: " +num2str(std(Results(:,i),'omitnan')) + newline + "N : " +num2str(length(find(isnan(Results(:,i))==0)))) 
end


function ExportCSV(hMotilityGui)
global Config;

Results = getappdata(hMotilityGui.fig,'CompResults');
Titles = getappdata(hMotilityGui.fig,'ResultTitles');

[pathname, filename, ~] = fileparts(fullfile(Config.Directory, Config.StackName));
savefile = fullfile(pathname, strcat(filename, '.xlsx'));
% check if files already exist. If so, then append or make a new file.
if isfile(savefile)
    olddata = readmatrix(savefile);
    Results = [olddata; Results];
end

Results = [Titles; Results];
writematrix(Results, savefile)


function Cancel(hMotilityGui)
close(hMotilityGui.fig);



