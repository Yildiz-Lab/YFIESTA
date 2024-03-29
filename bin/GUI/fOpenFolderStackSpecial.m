function output = fOpenFolderStackSpecial

c = 1.2*get(gcf,'Color');
ctext = 0.9*ones(1,3);
cbutton = [0.16,0.16,0.16];
ctextbox = [0.12,0.12,0.12];

hOpenSpecial.fig = figure('Units','normalized','WindowStyle','normal','DockControls','off','IntegerHandle','off','MenuBar','none','Name','Open Folder Stack Special',...
                      'NumberTitle','off','Position',[0.65 0.15 0.35 0.7],'HandleVisibility','callback','Tag','hOpenSpecial',...
                      'Visible','off','Resize','off','Color',c);
                  
fPlaceFig(hOpenSpecial.fig ,'special');

c = get(hOpenSpecial.fig ,'Color');
                  
hOpenSpecial.pMode = uibuttongroup('Parent',hOpenSpecial.fig,'Units','normalized','Position',[0.05 0.7 0.9 0.275],...
                                  'Title','Mode','Tag','tRange','FontSize',10,'SelectionChangeFcn',@ModeSelect,'BackgroundColor',c,'ForegroundColor',ctext);                  
                                  
% hOpenSpecial.rSeparateFiles = uicontrol('Parent',hOpenSpecial.pMode,'Units','normalized','Position',[0.05 0.675 0.9 0.25],'Enable','on','FontSize',12,...
%                                    'String','Separate files for each channel','Style','radiobutton','Tag','rSeparateFiles','HorizontalAlignment','left','BackgroundColor',c);
                                
hOpenSpecial.rSequentialSplitting = uicontrol('Parent',hOpenSpecial.pMode,'Units','normalized','Position',[0.05 0.375 0.9 0.25],'Enable','on','FontSize',12,...
                                    'String','Separate frames for each channel','Style','radiobutton','Tag','rSequentialSplitting','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);  
                                
hOpenSpecial.rSpatialSplitting = uicontrol('Parent',hOpenSpecial.pMode,'Units','normalized','Position',[0.05 0.075 0.9 0.25],'Enable','on','FontSize',12,...
                                    'String','Separate areas for each channel','Style','radiobutton','Tag','rSpatialSplitting','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);               
                  
% hOpenSpecial.pSeparateFiles = uipanel('Parent',hOpenSpecial.fig,'Units','normalized','Position',[0.05 0.275 0.9 0.4],'visible','on',...
%                                   'Title','Separate Files','Tag','pSeparateFiles','FontSize',10,'BackgroundColor',c);

% hOpenSpecial.pSeparate.tChannel(1) = uicontrol('Parent',hOpenSpecial.pSeparateFiles,'Units','normalized','Position',[0.05 0.78 0.075 0.15],'Enable','on','FontSize',12,...
%                                  'String','Ch1:','Style','text','Tag','tChannel1','HorizontalAlignment','left','BackgroundColor',c);                 
%                             
% hOpenSpecial.pSeparate.eChannel(1) = uicontrol('Parent',hOpenSpecial.pSeparateFiles,'Units','normalized','Position',[0.15 0.8 0.7 0.15],'Enable','inactive','FontSize',10,...
%                                 'String','','Style','edit','Tag','1','HorizontalAlignment','left','BackgroundColor','white',...
%                                 'ButtonDownFcn',@ChooseFile);         
%                             
% hOpenSpecial.pSeparate.bChannel(1) = uicontrol('Parent',hOpenSpecial.pSeparateFiles,'Units','normalized','Position',[0.875 0.8 0.1 0.15],'Enable','on','FontSize',10,...
%                                 'String','Load','Style','pushbutton','Tag','1','HorizontalAlignment','center',...
%                                 'Callback',@ChooseFile);  
%                             
% hOpenSpecial.pSeparate.tChannel(2) = uicontrol('Parent',hOpenSpecial.pSeparateFiles,'Units','normalized','Position',[0.05 0.53 0.075 0.15],'Enable','on','FontSize',12,...
%                                  'String','Ch2:','Style','text','Tag','tChannel2','HorizontalAlignment','left','BackgroundColor',c);                 
%                             
% hOpenSpecial.pSeparate.eChannel(2) = uicontrol('Parent',hOpenSpecial.pSeparateFiles,'Units','normalized','Position',[0.15 0.55 0.7 0.15],'Enable','inactive','FontSize',10,...
%                                 'String','','Style','edit','Tag','2','HorizontalAlignment','left','BackgroundColor','white',...
%                                 'ButtonDownFcn',@ChooseFile);         
%                             
% hOpenSpecial.pSeparate.bChannel(2) = uicontrol('Parent',hOpenSpecial.pSeparateFiles,'Units','normalized','Position',[0.875 0.55 0.1 0.15],'Enable','on','FontSize',10,...
%                                 'String','Load','Style','pushbutton','Tag','2','HorizontalAlignment','center',...
%                                 'Callback',@ChooseFile);  
%                             
% hOpenSpecial.pSeparate.tChannel(3) = uicontrol('Parent',hOpenSpecial.pSeparateFiles,'Units','normalized','Position',[0.05 0.28 0.075 0.15],'Enable','on','FontSize',12,...
%                                  'String','Ch3:','Style','text','Tag','tChannel3','HorizontalAlignment','left','BackgroundColor',c);                 
%                             
% hOpenSpecial.pSeparate.eChannel(3) = uicontrol('Parent',hOpenSpecial.pSeparateFiles,'Units','normalized','Position',[0.15 0.3 0.7 0.15],'Enable','inactive','FontSize',10,...
%                                 'String','','Style','edit','Tag','3','HorizontalAlignment','left','BackgroundColor','white',...
%                                 'ButtonDownFcn',@ChooseFile);         
%                             
% hOpenSpecial.pSeparate.bChannel(3) = uicontrol('Parent',hOpenSpecial.pSeparateFiles,'Units','normalized','Position',[0.875 0.3 0.1 0.15],'Enable','on','FontSize',10,...
%                                 'String','Load','Style','pushbutton','Tag','3','HorizontalAlignment','center',...
%                                 'Callback',@ChooseFile); 
%                             
% hOpenSpecial.pSeparate.tChannel(4) = uicontrol('Parent',hOpenSpecial.pSeparateFiles,'Units','normalized','Position',[0.05 0.03 0.075 0.15],'Enable','on','FontSize',12,...
%                                  'String','Ch4:','Style','text','Tag','tChannel4','HorizontalAlignment','left','BackgroundColor',c);                 
%                             
% hOpenSpecial.pSeparate.eChannel(4) = uicontrol('Parent',hOpenSpecial.pSeparateFiles,'Units','normalized','Position',[0.15 0.05 0.7 0.15],'Enable','inactive','FontSize',10,...
%                                 'String','','Style','edit','Tag','4','HorizontalAlignment','left','BackgroundColor','white',...
%                                 'ButtonDownFcn',@ChooseFile);         
%                             
% hOpenSpecial.pSeparate.bChannel(4) = uicontrol('Parent',hOpenSpecial.pSeparateFiles,'Units','normalized','Position',[0.875 0.05 0.1 0.15],'Enable','on','FontSize',10,...
%                                 'String','Load','Style','pushbutton','Tag','4','HorizontalAlignment','center',...
%                                 'Callback',@ChooseFile); 
     
hOpenSpecial.pSequentialSplitting = uipanel('Parent',hOpenSpecial.fig,'Units','normalized','Position',[0.05 0.275 0.9 0.4],'visible','on',...
                                     'Title','Sequential Frames','Tag','pSequentialSplitting','FontSize',10,'BackgroundColor',c,'ForegroundColor',ctext);
  
hOpenSpecial.pSequential.tNumberChannels = uicontrol('Parent',hOpenSpecial.pSequentialSplitting,'Units','normalized','Position',[0.05 0.8 0.5 0.15],'Enable','on','FontSize',12,...
                                   'String','Number of channels:','Style','text','Tag','tNumberChannels','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);
                               
hOpenSpecial.pSequential.mNumberChannels = uicontrol('Parent',hOpenSpecial.pSequentialSplitting,'Units','normalized','Position',[0.6 0.825 0.3 0.15],'Enable','on','FontSize',12,...
                                   'String',{'2','3','4'},'Style','popupmenu','Tag','mNumberChannels','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext,'Callback',@NumberSelect);
                               
hOpenSpecial.pSequential.rAlternateFrames = uicontrol('Parent',hOpenSpecial.pSequentialSplitting,'Units','normalized','Position',[0.05 0.65 0.9 0.15],'Value',1,'Enable','on','FontSize',12,...
                                   'Callback',@SequentialSelect,'String','Alternating frames','Style','radiobutton','Tag','rAlternateFrames','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);
                                 
hOpenSpecial.pSequential.rAlternateBlocks = uicontrol('Parent',hOpenSpecial.pSequentialSplitting,'Units','normalized','Position',[0.05 0.45 0.9 0.15],'Enable','on','FontSize',12,...
                                   'Callback',@SequentialSelect,'String','Alternating blocks','Style','radiobutton','Tag','rAlternateBlocks','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);
                               
hOpenSpecial.pSequential.rSingleBlocks = uicontrol('Parent',hOpenSpecial.pSequentialSplitting,'Units','normalized','Position',[0.05 0.25 0.9 0.15],'Enable','on','FontSize',12,...
                                   'Callback',@SequentialSelect,'String','Single blocks','Style','radiobutton','Tag','rSingleBlocks','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);
                               
hOpenSpecial.pSequential.tCh(1) = uicontrol('Parent',hOpenSpecial.pSequentialSplitting,'Units','normalized','Position',[0.05 0.02 0.09 0.15],'Enable','on','FontSize',12,...
                                 'String','Ch1:','Style','text','Tag','tCh1','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);                 
                            
hOpenSpecial.pSequential.eCh(1) = uicontrol('Parent',hOpenSpecial.pSequentialSplitting,'Units','normalized','Position',[0.15 0.05 0.1 0.15],'Enable','off','FontSize',10,...
                                'String','','Style','edit','Tag','eCh1','HorizontalAlignment','center','BackgroundColor',ctextbox,'ForegroundColor',ctext);         
                                                
hOpenSpecial.pSequential.tCh(2) = uicontrol('Parent',hOpenSpecial.pSequentialSplitting,'Units','normalized','Position',[0.275 0.02 0.09 0.15],'Enable','on','FontSize',12,...
                                 'String','Ch2:','Style','text','Tag','tCh2','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);                 
                            
hOpenSpecial.pSequential.eCh(2) = uicontrol('Parent',hOpenSpecial.pSequentialSplitting,'Units','normalized','Position',[0.375 0.05 0.1 0.15],'Enable','off','FontSize',10,...
                                'String','','Style','edit','Tag','eCh2','HorizontalAlignment','center','BackgroundColor',ctextbox,'ForegroundColor',ctext);         
                                                    
hOpenSpecial.pSequential.tCh(3) = uicontrol('Parent',hOpenSpecial.pSequentialSplitting,'Units','normalized','Position',[0.5 0.02 0.09 0.15],'Enable','on','FontSize',12,...
                                 'String','Ch3:','Style','text','Tag','tCh3','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);                 
                            
hOpenSpecial.pSequential.eCh(3) = uicontrol('Parent',hOpenSpecial.pSequentialSplitting,'Units','normalized','Position',[0.6 0.05 0.1 0.15],'Enable','off','FontSize',10,...
                                'String','','Style','edit','Tag','eCh3','HorizontalAlignment','center','BackgroundColor',ctextbox,'ForegroundColor',ctext);         
 
hOpenSpecial.pSequential.tCh(4) = uicontrol('Parent',hOpenSpecial.pSequentialSplitting,'Units','normalized','Position',[0.725 0.02 0.09 0.15],'Enable','on','FontSize',12,...
                                 'String','Ch4:','Style','text','Tag','tCh4','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);                 
                            
hOpenSpecial.pSequential.eCh(4) = uicontrol('Parent',hOpenSpecial.pSequentialSplitting,'Units','normalized','Position',[0.825 0.05 0.1 0.15],'Enable','off','FontSize',10,...
                                'String','','Style','edit','Tag','eCh4','HorizontalAlignment','center','BackgroundColor',ctextbox,'ForegroundColor',ctext);         
                                                        
hOpenSpecial.pSpatialSplitting = uibuttongroup('Parent',hOpenSpecial.fig,'Units','normalized','Position',[0.05 0.275 0.9 0.4],'visible','off',...
                                 'Title','Separate areas','Tag','pSpatialSplitting','FontSize',10,'BackgroundColor',c,'ForegroundColor',ctext);   

hOpenSpecial.pSpatial.rTwoHorizontal = uicontrol('Parent',hOpenSpecial.pSpatialSplitting,'Units','normalized','Position',[0.05 0.8 0.9 0.15],'Enable','on','FontSize',12,...
                                        'String','Two channels, split horizontal','Style','radiobutton','Tag','rTwoHorizontal','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);
                                    
                                    
hOpenSpecial.pSpatial.rTwoVertical = uicontrol('Parent',hOpenSpecial.pSpatialSplitting,'Units','normalized','Position',[0.05 0.625 0.9 0.15],'Enable','on','FontSize',12,...
                                        'String','Two channels, spilt vertical','Style','radiobutton','Tag','rTwoVertical','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);                                    
                                    
% JS Edit 2022/07/27 for custom split since the gemini position might not
% be perfectly centrally aligned. This should take off the translation
% error

hOpenSpecial.pSpatial.rTwoHorizontalCustom1text = uicontrol('Parent',hOpenSpecial.pSpatialSplitting,'Units','normalized','Position',[0.55 0.715 0.1 0.15],'Enable','on','FontSize',12,...
                                        'String','center','Style','text','Tag','rTwoHorizontalCustom','BackgroundColor',c,'ForegroundColor',ctext);

hOpenSpecial.pSpatial.rTwoHorizontalCustom1 = uicontrol('Parent',hOpenSpecial.pSpatialSplitting,'Units','normalized','Position',[0.65 0.715 0.1 0.15],'Enable','on','FontSize',12,...
                                        'String','525','Style','edit','Tag','rTwoHorizontalCustom1','BackgroundColor',ctextbox,'ForegroundColor',ctext);
     
hOpenSpecial.pSpatial.rTwoHorizontalCustom2text = uicontrol('Parent',hOpenSpecial.pSpatialSplitting,'Units','normalized','Position',[0.77 0.715 0.1 0.15],'Enable','on','FontSize',12,...
                                        'String','width','Style','text','Tag','rTwoHorizontalCustom','BackgroundColor',c,'ForegroundColor',ctext);

hOpenSpecial.pSpatial.rTwoHorizontalCustom2 = uicontrol('Parent',hOpenSpecial.pSpatialSplitting,'Units','normalized','Position',[0.87 0.715 0.1 0.15],'Enable','on','FontSize',12,...
                                        'String','498','Style','edit','Tag','rTwoHorizontalCustom2','BackgroundColor',ctextbox,'ForegroundColor',ctext);
% End of JS Edit

                                    
hOpenSpecial.pSpatial.rFourSymmetric = uicontrol('Parent',hOpenSpecial.pSpatialSplitting,'Units','normalized','Position',[0.05 0.45 0.9 0.15],'Enable','on','FontSize',12,...
                                        'String','Four channels, 2x2 symmetric','Style','radiobutton','Tag','rFourSymmetric','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);                              
                                 
hOpenSpecial.pSpatial.rFourHorizontal = uicontrol('Parent',hOpenSpecial.pSpatialSplitting,'Units','normalized','Position',[0.05 0.275 0.9 0.15],'Enable','on','FontSize',12,...
                                        'String','Four channels, split horizontal','Style','radiobutton','Tag','rFourHorizontal','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext); 
    
hOpenSpecial.pSpatial.rFourHorizontal = uicontrol('Parent',hOpenSpecial.pSpatialSplitting,'Units','normalized','Position',[0.05 0.1 0.9 0.15],'Enable','on','FontSize',12,...
                                        'String','Four channels, split vertical','Style','radiobutton','Tag','rFourVertical','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext); 

hOpenSpecial.pPostProcessing = uipanel('Parent',hOpenSpecial.fig,'Units','normalized','Position',[0.05 0.05 0.6 0.2],...
                                     'Title','Optional Postprocessing','Tag','pPostProcessing','FontSize',10,'BackgroundColor',c,'ForegroundColor',ctext);          
                                 
hOpenSpecial.pPost.cParallax = uicontrol('Parent',hOpenSpecial.pPostProcessing,'Units','normalized','Position',[0.05 0.725 0.9 0.25],'Enable','on','FontSize',12,...
                                        'Callback',@PostSelect,'String','Parallax 3D Tracking','Style','checkbox','Tag','cParallax','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);     
                                    
hOpenSpecial.pPost.cPolTIRF = uicontrol('Parent',hOpenSpecial.pPostProcessing,'Units','normalized','Position',[0.05 0.4 0.9 0.25],'Enable','on','FontSize',12,...
                                        'Callback',@PostSelect,'String','PolTIRF Orientation Tracking','Style','checkbox','Tag','cPolTIRF','HorizontalAlignment','left','BackgroundColor',c,'ForegroundColor',ctext);  

                                    % JS Edit 2022/09/09 for summing frames
                                    % (antiblock)
hOpenSpecial.pPost.rAntiblockCh1text = uicontrol('Parent',hOpenSpecial.pPostProcessing,'Units','normalized','Position',[0.02 0.1 0.15 0.25],'Enable','on','FontSize',10,...
                                        'String','Ch1','Style','text','Tag','rAntiblockCh1text','BackgroundColor',c,'ForegroundColor',ctext);
                                    
hOpenSpecial.pPost.rAntiblockCh1 = uicontrol('Parent',hOpenSpecial.pPostProcessing,'Units','normalized','Position',[0.15 0.1 0.1 0.25],'Enable','on','FontSize',12,...
                                        'String','1','Style','edit','Tag','rAntiblockCh1','BackgroundColor',ctextbox,'ForegroundColor',ctext);

hOpenSpecial.pPost.rAntiblockCh2text = uicontrol('Parent',hOpenSpecial.pPostProcessing,'Units','normalized','Position',[0.25 0.1 0.15 0.25],'Enable','on','FontSize',10,...
                                        'String','Ch2','Style','text','Tag','rAntiblockCh2text','BackgroundColor',c,'ForegroundColor',ctext);
                                    
hOpenSpecial.pPost.rAntiblockCh2 = uicontrol('Parent',hOpenSpecial.pPostProcessing,'Units','normalized','Position',[0.4 0.1 0.1 0.25],'Enable','on','FontSize',12,...
                                        'String','1','Style','edit','Tag','rAntiblockCh2','BackgroundColor',ctextbox,'ForegroundColor',ctext);                                    

hOpenSpecial.pPost.rAntiblockCh3text = uicontrol('Parent',hOpenSpecial.pPostProcessing,'Units','normalized','Position',[0.5 0.1 0.15 0.25],'Enable','on','FontSize',10,...
                                        'String','Ch3','Style','text','Tag','rAntiblockCh3text','BackgroundColor',c,'ForegroundColor',ctext);
                                    
hOpenSpecial.pPost.rAntiblockCh3 = uicontrol('Parent',hOpenSpecial.pPostProcessing,'Units','normalized','Position',[0.65 0.1 0.1 0.25],'Enable','on','FontSize',12,...
                                        'String','1','Style','edit','Tag','rAntiblockCh3','BackgroundColor',ctextbox,'ForegroundColor',ctext);  
                                    
hOpenSpecial.pPost.rAntiblockCh4text = uicontrol('Parent',hOpenSpecial.pPostProcessing,'Units','normalized','Position',[0.75 0.1 0.15 0.25],'Enable','on','FontSize',10,...
                                        'String','Ch4','Style','text','Tag','rAntiblockCh4text','BackgroundColor',c,'ForegroundColor',ctext);
                                    
hOpenSpecial.pPost.rAntiblockCh4 = uicontrol('Parent',hOpenSpecial.pPostProcessing,'Units','normalized','Position',[0.9 0.1 0.1 0.25],'Enable','on','FontSize',12,...
                                        'String','1','Style','edit','Tag','rAntiblockCh4','BackgroundColor',ctextbox,'ForegroundColor',ctext);  

                                    % End of JS Edit 2022/09/09 
                                    
hOpenSpecial.bLoadStack = uicontrol('Parent',hOpenSpecial.fig,'Units','normalized','Position',[0.7 0.05 0.25 0.15],'Enable','on','FontSize',14,...
                                'String','Load Stack','Style','pushbutton','Tag','bLoadStack','HorizontalAlignment','center','BackgroundColor',cbutton,'ForegroundColor',ctext,...
                                'Callback',@LoadStack);  

%set(hOpenSpecial.rSeparateFiles,'UserData',hOpenSpecial.pSeparateFiles);
set(hOpenSpecial.rSequentialSplitting,'UserData',hOpenSpecial.pSequentialSplitting);
set(hOpenSpecial.rSpatialSplitting,'UserData',hOpenSpecial.pSpatialSplitting);

setappdata(0,'hOpenSpecial',hOpenSpecial);
setappdata(hOpenSpecial.fig,'output',[]);
uiwait(hOpenSpecial.fig)
try
    output = getappdata(hOpenSpecial.fig,'output');
    close(hOpenSpecial.fig);
catch 
    output =[];
end

function ModeSelect(~,o)
set(get(o.NewValue,'UserData'),'visible','on');
set(get(o.OldValue,'UserData'),'visible','off');

function SequentialSelect(~,o)
hOpenSpecial = getappdata(0,'hOpenSpecial');
set(findobj('Parent',get(o.Source,'Parent'),'-and','Style','radiobutton'),'Value',0);
set(o.Source,'Value',1);
n = get(hOpenSpecial.pSequential.mNumberChannels,'Value');
h = findobj('Parent',get(o.Source,'Parent'),'-and','Style','edit');
if strcmp(get(o.Source,'Tag'),'rAlternateFrames')
    set(h,'enable','off','String','');
elseif strcmp(get(o.Source,'Tag'),'rAlternateBlocks')
    set(h((4-n):4),'enable','on','String','10');
    set(h(1:(3-n)),'enable','off','String','');
else
    set(h((5-n):4),'enable','on','String','1');
    set(h((4-n)),'enable','inactive','String','rest');
    set(h(1:(3-n)),'enable','off','String','');
end

function PostSelect(h,o)
if h.Value == 1
    set(findobj('Parent',get(o.Source,'Parent'),'-and','Style','checkbox'),'Value',0);
    set(o.Source,'Value',1);   
end

function NumberSelect(h,o)
n = h.Value;
h = findobj('Parent',get(o.Source,'Parent'),'-and','Style','radiobutton');
c = get(h,'Value');
h = findobj('Parent',get(o.Source,'Parent'),'-and','Style','edit');
if c{3} == 1
    set(h,'enable','off','String','');
elseif c{2} == 1
    set(h((4-n):4),'enable','on','String','10');
    set(h(1:(3-n)),'enable','off','String','');
else
    set(h((5-n):4),'enable','on','String','1');
    set(h((4-n)),'enable','inactive','String','rest');
    set(h(1:(3-n)),'enable','off','String','');
end

% function ChooseFile(~,o)
% global FiestaDir;
% hOpenSpecial = getappdata(0,'hOpenSpecial');
% [FileName,PathName] = uigetfile({'*.stk;*.nd;*.nd2;*.dv;*.zvi;*.tif;*.tiff','Image Stacks (*.stk,*.nd,*.nd2,*.dv,*.zvi,*.tif,*.tiff)'},'Select the Stack',FiestaDir.Stack,'MultiSelect','on'); %open dialog for *.stk files
% if PathName~=0
%     FiestaDir.Stack=PathName;
%     n = str2double(get(o.Source,'Tag'));
%     if ~iscell(FileName)
%         FileName ={FileName};
%     end
%     if length(FileName)+n>5
%         e = 4-n+1;
%     else
%         e = length(FileName);
%     end
%     for m = 1:e
%         set(hOpenSpecial.pSeparate.eChannel(m+n-1),'String',FileName{m},'UserData',PathName);
%     end
% end

% function ChooseFile(~,o) %actually folder
% global FiestaDir;
% hOpenSpecial = getappdata(0,'hOpenSpecial');
% %[FileName,PathName] = uigetfile({'*.stk;*.nd;*.nd2;*.dv;*.zvi;*.tif;*.tiff','Image Stacks (*.stk,*.nd,*.nd2,*.dv,*.zvi,*.tif,*.tiff)'},'Select the Stack',FiestaDir.Stack,'MultiSelect','on'); %open dialog for *.stk files
% PathName = uigetdir('Image Folder or Superfolder');
% if PathName~=0
%     FiestaDir.Stack=PathName;
%     n = str2double(get(o.Source,'Tag'));
% %     if ~iscell(FileName)
% %         FileName ={FileName};
% %     end
%     if length(PathName)+n>5
%         e = 4-n+1;
%     else 
%         e = length(PathName);
%     end
%     for m = 1:e
%         set(hOpenSpecial.pSeparate.eChannel(m+n-1),'String',PathName{m},'UserData',PathName);
%     end
% end

function LoadStack(~,~)
global FiestaDir;
hOpenSpecial = getappdata(0,'hOpenSpecial');
output.Mode = get(get(hOpenSpecial.pMode,'SelectedObject'),'Tag');
if strcmp(output.Mode,'rSeparateFiles')
    output.Data = get(hOpenSpecial.pSeparate.eChannel,{'String','UserData'});
    for n = size(output.Data,1):-1:1
        if isempty(output.Data{n,1})
            output.Data(n,:) = [];
        end
    end
elseif strcmp(output.Mode,'rSequentialSplitting')
    %[FileName,PathName] = uigetfile({'*.stk;*.nd;*.nd2;*.dv;*.zvi;*.tif;*.tiff','Image Stacks (*.stk,*.nd,*.nd2,*.dv,*.zvi,*.tif,*.tiff)'},'Select the Stack',FiestaDir.Stack,'MultiSelect','on'); %open dialog for *.stk files
    PathName = uigetdir('Image Folder or Superfolder');
    FileName = '';
    if PathName~=0
        FiestaDir.Stack=PathName;
        output.Data = {FileName,PathName,get(hOpenSpecial.pSequential.mNumberChannels,'Value')+1,...
                   get(findobj('Parent',hOpenSpecial.pSequentialSplitting,'-and','Style','radiobutton'),'Value'),get(hOpenSpecial.pSequential.eCh,'String')};
    end
else
    %[FileName,PathName] = uigetfile({'*.stk;*.nd;*.nd2;*.dv;*.zvi;*.tif;*.tiff','Image Stacks (*.stk,*.nd,*.nd2,*.dv,*.zvi,*.tif,*.tiff)'},'Select the Stack',FiestaDir.Stack,'MultiSelect','on'); %open dialog for *.stk files
    PathName = uigetdir('Image Folder or Superfolder');
    FileName = '';
    if PathName~=0
        FiestaDir.Stack=PathName;
        % JS Edit 2022/07/27 to load in custom string for horizontal
        % splitting
        c=get(findobj('Parent',hOpenSpecial.pSpatialSplitting,'Tag','rTwoHorizontalCustom1'),'String');
        w=get(findobj('Parent',hOpenSpecial.pSpatialSplitting,'Tag','rTwoHorizontalCustom2'),'String');
        output.Data = {FileName,PathName,get(findobj('Parent',hOpenSpecial.pSpatialSplitting,'-and','Style','radiobutton'),'Value'),[str2double(c),str2double(w)]};
    end
end
ch = get(hOpenSpecial.pPostProcessing,'Children');
select = cell2mat(get(ch,'Value'));
output.Optional = get(ch(select==1),'Tag');
if ~isempty(output.Optional)
    output.Optional(1) = [];    
end

%JS Edit 2022/09/09 for antiblock
ch1 = get(findobj('Parent',hOpenSpecial.pPostProcessing,'Tag','rAntiblockCh1'),'String');
ch2 = get(findobj('Parent',hOpenSpecial.pPostProcessing,'Tag','rAntiblockCh2'),'String');
ch3 = get(findobj('Parent',hOpenSpecial.pPostProcessing,'Tag','rAntiblockCh3'),'String');
ch4 = get(findobj('Parent',hOpenSpecial.pPostProcessing,'Tag','rAntiblockCh4'),'String');
output.PostSumFrames = {str2num(ch1), str2num(ch2), str2num(ch3), str2num(ch4)};
%End of JS Edit 2022/09/09 for antiblock

setappdata(hOpenSpecial.fig,'output',output);
uiresume(gcbf);