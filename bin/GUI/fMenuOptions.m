function fMenuOptions(func,varargin)
switch func
    case 'LoadConfig'
        LoadConfig(varargin{1});
    case 'SaveConfig'
        SaveConfig;
    case 'SetDefaultConfig'
        SetDefaultConfig;
     case 'SaveCorrections'
        SaveCorrections(varargin{1});
    case 'LoadCorrections'
        LoadCorrections(varargin{1});
    case 'SetMaskCorrections'
        SetMaskCorrections(varargin{1})
    case 'ApplyMultiChannelCorrections'
        ApplyMultiChannelCorrections(varargin{1});
end

function LoadConfig(hMainGui)
global Config;
[FileName, PathName] = uigetfile({'*.mat','FIESTA Config(*.mat)'},'Load FIESTA Config',fShared('GetLoadDir'));
if FileName~=0
    fShared('SetLoadDir',PathName);
    tempConfig=fLoad([PathName FileName],'Config');
    Config.ConnectMol=tempConfig.ConnectMol;
    Config.ConnectFil=tempConfig.ConnectFil;    
    Config.Threshold=tempConfig.Threshold;
    Config.RefPoint=tempConfig.RefPoint;
    Config.OnlyTrack=tempConfig.OnlyTrack;
    Config.BorderMargin=tempConfig.BorderMargin;
    Config.Model=tempConfig.Model;
    Config.MaxFunc=tempConfig.MaxFunc;
    Config.DynamicFil=tempConfig.DynamicFil;
    Config.ReduceFitBox=tempConfig.ReduceFitBox;
    Config.FilFocus=tempConfig.FilFocus;
end
fShow('Image',hMainGui);

function SaveConfig
global Config; %#ok<NUSED>
[FileName, PathName] = uiputfile({'*.mat','MAT-File(*.mat)'},'Save FIESTA Config',fShared('GetSaveDir'));
if FileName~=0
    fShared('SetSaveDir',PathName);
    file = [PathName FileName];
    if isempty(findstr('.mat',file))
        file = [file '.mat'];
    end
    save(file,'Config');
end

function SetDefaultConfig
global Config;
global DirCurrent
button = fQuestDlg('Overwrite the default configuration?','FIESTA Warning',{'Overwrite','Cancel'},'Cancel');
if strcmp(button,'Overwrite')==1
    if isdeployed
        if ismac
            file_id = fopen('~/Library/Fiesta/fiesta.ini','w');
        elseif ispc
            file_id = fopen([winqueryreg('HKEY_CURRENT_USER','Software\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders','Local AppData') '\Fiesta\fiesta.ini'],'w');
        end
    else
        file_id = fopen([DirCurrent 'fiesta.ini'],'w');
    end
    fwrite(file_id,jsonencode(Config));
    fclose(file_id);
end

function LoadCorrections(hMainGui)
global Stack;
fRightPanel('CheckReference',hMainGui);
[FileName, PathName] = uigetfile({'*.mat','FIESTA Transformation (*.mat)'},'Load FIESTA Reference Transformations',fShared('GetLoadDir'));
if FileName~=0
    fShared('SetLoadDir',PathName);    
    Drift=fLoad([PathName FileName],'Drift');
    if ~isempty(Drift)
        if ~iscell(Drift)
            fMsgDlg('References not compatible with this FIESTA version','error');
            return;
        end
        if numel(Stack)>numel(Drift)
            Drift{numel(Stack)} = [];
        end
        setappdata(hMainGui.fig,'Drift',Drift);
    end
    fShared('UpdateMenu',hMainGui);
end
setappdata(0,'hMainGui',hMainGui);

% JS Edit 2022/09/03
% Putting the ability to apply corrections to Ch2 Molecules
function SetMaskCorrections(hMainGui)

LoadDir = fShared('GetLoadDir');  % In FIESTA
[baseName, folder] = uigetfile({'*.mat','FIESTA Data(*.mat)'},'Load FIESTA Objects for Channel 1',LoadDir,'MultiSelect','off');
ch1file = fullfile(folder, baseName);
if strcmp(folder, LoadDir) == 0
    LoadDir = folder;
end
[baseName, folder] = uigetfile({'*.mat','FIESTA Data(*.mat)'},'Load FIESTA Objects for Channel 2',LoadDir,'MultiSelect','off');
ch2file = fullfile(folder, baseName);

F = fCalcAlignCorrection(ch1file, ch2file);
setappdata(hMainGui.fig,'MaskCorrect',F);
setappdata(hMainGui.fig,'MaskCorrectFiles',[ch1file, ch2file]);
setappdata(hMainGui.fig,'MaskCorrectionApplied',0);

set(hMainGui.Menu.mApplyMultiCorrections,'Enable','On');

setappdata(0,'hMainGui',hMainGui);

function ApplyMultiChannelCorrections(hMainGui)
global Stack
global Molecule

F = getappdata(hMainGui.fig,'MaskCorrect');
Fx = F{1}; Fy = F{2};
toggle = getappdata(hMainGui.fig,'MaskCorrectionApplied');

if toggle %we already corrected, go back
    % maybe a dialog box would be helpful
    answer = questdlg('Transformation already applied, apply inverse transformation?');
    if length(answer) < 3 % No
        return
    elseif length(answer) < 4 %Yes
        pm = -1;
    else %Cancel
        return
    end
else %not corrected so do normal transformation
    pm = 1;
end

for i = 1:length(Molecule)
    if Molecule(i).Channel == 2
        % for some reason it needs to be a double array
        Molecule(i).Results(:,3) = Molecule(i).Results(:,3) + pm*Fx(double(Molecule(i).Results(:,3:4)));
        Molecule(i).Results(:,4) = Molecule(i).Results(:,4) + pm*Fy(double(Molecule(i).Results(:,3:4)));
        
    end
    
    Molecule(i).Results(:,6) = sqrt( (Molecule(i).Results(:,3) - Molecule(i).Results(1,3)).^2 + (Molecule(i).Results(:,4) - Molecule(i).Results(1,4)).^2 );
    
end

toggle = mod(toggle+1,2);
setappdata(hMainGui.fig,'MaskCorrectionApplied',toggle)

% the global variable enscribes that this change happens imediately in the
% Results. But we still want it to plot again so use built in function to
% show tracks

fShow('Tracks')
% End of JS Edit


function SaveCorrections(hMainGui)
Drift=getappdata(hMainGui.fig,'Drift'); %#ok<NASGU>
[FileName, PathName] = uiputfile({'*.mat','MAT-files (*.mat)'},'Save FIESTA Reference Transformations',fShared('GetSaveDir'));
if FileName~=0
    fShared('SetSaveDir',PathName);
    file = strtok(FileName,'.');
    save([PathName file '.mat'],'Drift');  
end

