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
    case 'SetLinearMaskCorrections'
        SetLinearMaskCorrections(varargin{1});
    case 'ApplyMultiChannelCorrections'
        ApplyMultiChannelCorrections(varargin{1});
    case 'SaveMaskCorrections'
        SaveMaskCorrections(varargin{1});
    case 'LoadMaskCorrections'
        LoadMaskCorrections(varargin{1});
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
% DC 260721 rewritten to pick which channels the corrections are applied to
function SetMaskCorrections(hMainGui)
global Molecule;

% If a previous correction is still applied, revert it first so we start
% from clean/raw positions before computing the new one
prevApplied = getappdata(hMainGui.fig,'MaskCorrectionApplied');
if isequal(prevApplied,1)
    prevDeltas = getappdata(hMainGui.fig,'MaskCorrectionDeltas');
    prevChannels = getappdata(hMainGui.fig,'MaskCorrectChannels');
    if ~isempty(prevDeltas) && ~isempty(prevChannels)
        prevAdjChannel = prevChannels(2);
        for i = 1:length(Molecule)
            if Molecule(i).Channel == prevAdjChannel && i <= length(prevDeltas) && ~isempty(prevDeltas{i})
                Molecule(i).Results(:,3) = Molecule(i).Results(:,3) + prevDeltas{i}(:,1);
                Molecule(i).Results(:,4) = Molecule(i).Results(:,4) + prevDeltas{i}(:,2);
                Molecule(i).Results(:,6) = sqrt( (Molecule(i).Results(:,3) - Molecule(i).Results(1,3)).^2 + (Molecule(i).Results(:,4) - Molecule(i).Results(1,4)).^2 );
            end
        end
    end
    setappdata(hMainGui.fig,'MaskCorrectionApplied',0);
end

refChannel = str2double(fInputDlg('Enter reference channel number (stays fixed):','1'));
if isnan(refChannel)
    return;
end
adjChannel = str2double(fInputDlg('Enter channel number to adjust:','2'));
if isnan(adjChannel)
    return;
end
too_large_nm = str2double(fInputDlg('Max matching distance between channels (nm). Increase if objects further apart than this are being missed:','100'));
if isnan(too_large_nm) || too_large_nm <= 0
    return;
end

LoadDir = fShared('GetLoadDir');
[baseName, folder] = uigetfile({'*.mat','FIESTA Data(*.mat)'}, ...
    sprintf('Load FIESTA Objects for Channel %d (reference)',refChannel),LoadDir,'MultiSelect','off');
if isequal(baseName,0)
    return;
end
ch1file = fullfile(folder, baseName);
fShared('SetLoadDir',folder);

[baseName, folder] = uigetfile({'*.mat','FIESTA Data(*.mat)'}, ...
    sprintf('Load FIESTA Objects for Channel %d (to adjust)',adjChannel),folder,'MultiSelect','off');
if isequal(baseName,0)
    return;
end
ch2file = fullfile(folder, baseName);
fShared('SetLoadDir',folder);

F = fCalcAlignCorrection(ch1file, ch2file, refChannel, adjChannel, too_large_nm);
Fx = F{1}; Fy = F{2};

deltas = cell(1,length(Molecule));
for i = 1:length(Molecule)
    if Molecule(i).Channel == adjChannel
        rawXY = double(Molecule(i).Results(:,3:4));
        dx = Fx(rawXY);
        dy = Fy(rawXY);
        deltas{i} = [dx, dy];
        Molecule(i).Results(:,3) = Molecule(i).Results(:,3) - dx;
        Molecule(i).Results(:,4) = Molecule(i).Results(:,4) - dy;
        Molecule(i).Results(:,6) = sqrt( (Molecule(i).Results(:,3) - Molecule(i).Results(1,3)).^2 + (Molecule(i).Results(:,4) - Molecule(i).Results(1,4)).^2 );
    end
end

setappdata(hMainGui.fig,'MaskCorrect',F);
setappdata(hMainGui.fig,'MaskCorrectFiles',{ch1file, ch2file});
setappdata(hMainGui.fig,'MaskCorrectChannels',[refChannel, adjChannel]);
setappdata(hMainGui.fig,'MaskCorrectionDeltas',deltas);
setappdata(hMainGui.fig,'MaskCorrectionApplied',1);

set(hMainGui.Menu.mApplyMultiCorrections,'Enable','On');
set(hMainGui.Menu.mSaveMaskCorrections,'Enable','On');

setappdata(0,'hMainGui',hMainGui);
fShow('Tracks');


function SetLinearMaskCorrections(hMainGui)
[dx,dy] = offset_by_raw_means_v2(); %get straight correction

% Now get limits of current view for application
interpolation_factor = 200;
xb = hMainGui.MidPanel.aView.XLim(1):interpolation_factor:hMainGui.MidPanel.aView.XLim(2);
yb = hMainGui.MidPanel.aView.YLim(1):interpolation_factor:hMainGui.MidPanel.aView.YLim(2);
[x,y] = meshgrid(xb, yb);
x = reshape(x,[],1); y = reshape(y,[],1);
Fx = scatteredInterpolant(x, y, double(dx)*ones(numel(x),1), 'natural'); %gotta make doubles for interpolant
Fy = scatteredInterpolant(x, y, double(dy)*ones(numel(y),1), 'natural');
F = {Fx, Fy};

setappdata(hMainGui.fig,'MaskCorrect',F);
setappdata(hMainGui.fig,'MaskCorrectionApplied',0);

set(hMainGui.Menu.mApplyMultiCorrections,'Enable','On');

setappdata(0,'hMainGui',hMainGui);

function ApplyMultiChannelCorrections(hMainGui)
global Molecule;

deltas = getappdata(hMainGui.fig,'MaskCorrectionDeltas');
channels = getappdata(hMainGui.fig,'MaskCorrectChannels');
adjChannel = channels(2);
toggle = getappdata(hMainGui.fig,'MaskCorrectionApplied');

if toggle
    sign_ = -1;   % currently on -> subtract the cached delta back off
else
    sign_ = 1;    % currently off -> re-add the cached delta
end

for i = 1:length(Molecule)
    if Molecule(i).Channel == adjChannel && i <= length(deltas) && ~isempty(deltas{i})
        Molecule(i).Results(:,3) = Molecule(i).Results(:,3) + sign_*deltas{i}(:,1);
        Molecule(i).Results(:,4) = Molecule(i).Results(:,4) + sign_*deltas{i}(:,2);
        Molecule(i).Results(:,6) = sqrt( (Molecule(i).Results(:,3) - Molecule(i).Results(1,3)).^2 + (Molecule(i).Results(:,4) - Molecule(i).Results(1,4)).^2 );
    end
end

setappdata(hMainGui.fig,'MaskCorrectionApplied', ~toggle);

fShow('Tracks')
% end of DC edit
% End of JS Edit

% DC260721 added saving / loading multichannel corrections
function SaveMaskCorrections(hMainGui)
F = getappdata(hMainGui.fig,'MaskCorrect');
MaskCorrectFiles = getappdata(hMainGui.fig,'MaskCorrectFiles'); %#ok<NASGU>
MaskCorrectChannels = getappdata(hMainGui.fig,'MaskCorrectChannels'); %#ok<NASGU>

if isempty(F)
    fMsgDlg('No multichannel correction has been set yet.','error');
    return;
end

[FileName, PathName] = uiputfile({'*.mat','FIESTA Multichannel Correction (*.mat)'},'Save Multichannel Correction',fShared('GetSaveDir'));
if FileName~=0
    fShared('SetSaveDir',PathName);
    file = strtok(FileName,'.');
    save([PathName file '.mat'],'F','MaskCorrectFiles','MaskCorrectChannels');
end


function LoadMaskCorrections(hMainGui)
global Molecule;

[FileName, PathName] = uigetfile({'*.mat','FIESTA Multichannel Correction (*.mat)'},'Load Multichannel Correction',fShared('GetLoadDir'));
if isequal(FileName,0)
    return;
end
fShared('SetLoadDir',PathName);

loaded = load(fullfile(PathName,FileName),'F','MaskCorrectChannels');
if ~isfield(loaded,'F') || ~isfield(loaded,'MaskCorrectChannels')
    fMsgDlg('Selected file does not contain a valid multichannel correction.','error');
    return;
end
F = loaded.F;
channels = loaded.MaskCorrectChannels;
adjChannel = channels(2);

% revert any correction currently applied, same guard as SetMaskCorrections
prevApplied = getappdata(hMainGui.fig,'MaskCorrectionApplied');
if isequal(prevApplied,1)
    prevDeltas = getappdata(hMainGui.fig,'MaskCorrectionDeltas');
    prevChannels = getappdata(hMainGui.fig,'MaskCorrectChannels');
    if ~isempty(prevDeltas) && ~isempty(prevChannels)
        prevAdjChannel = prevChannels(2);
        for i = 1:length(Molecule)
            if Molecule(i).Channel == prevAdjChannel && i <= length(prevDeltas) && ~isempty(prevDeltas{i})
                Molecule(i).Results(:,3) = Molecule(i).Results(:,3) + prevDeltas{i}(:,1);
                Molecule(i).Results(:,4) = Molecule(i).Results(:,4) + prevDeltas{i}(:,2);
                Molecule(i).Results(:,6) = sqrt( (Molecule(i).Results(:,3) - Molecule(i).Results(1,3)).^2 + (Molecule(i).Results(:,4) - Molecule(i).Results(1,4)).^2 );
            end
        end
    end
    setappdata(hMainGui.fig,'MaskCorrectionApplied',0);
end

Fx = F{1}; Fy = F{2};

% compute all deltas first, validate, then apply (same pattern as SetMaskCorrections)
deltas = cell(1,length(Molecule));
badIdx = [];
for i = 1:length(Molecule)
    if Molecule(i).Channel == adjChannel
        rawXY = double(Molecule(i).Results(:,3:4));
        dx = Fx(rawXY);
        dy = Fy(rawXY);
        deltas{i} = [dx, dy];
        if any(isnan(dx)) || any(isnan(dy))
            badIdx(end+1) = i; %#ok<AGROW>
        end
    end
end

if ~isempty(badIdx)
    fMsgDlg({sprintf('%d molecule(s) fall outside the loaded correction''s matched domain.',numel(badIdx)),...
             'No correction was applied.'},'error');
    return;
end

for i = 1:length(Molecule)
    if Molecule(i).Channel == adjChannel
        Molecule(i).Results(:,3) = Molecule(i).Results(:,3) - deltas{i}(:,1);
        Molecule(i).Results(:,4) = Molecule(i).Results(:,4) - deltas{i}(:,2);
        Molecule(i).Results(:,6) = sqrt( (Molecule(i).Results(:,3) - Molecule(i).Results(1,3)).^2 + (Molecule(i).Results(:,4) - Molecule(i).Results(1,4)).^2 );
    end
end

setappdata(hMainGui.fig,'MaskCorrect',F);
setappdata(hMainGui.fig,'MaskCorrectChannels',channels);
setappdata(hMainGui.fig,'MaskCorrectionDeltas',deltas);
setappdata(hMainGui.fig,'MaskCorrectionApplied',1);

set(hMainGui.Menu.mApplyMultiCorrections,'Enable','On');
set(hMainGui.Menu.mSaveMaskCorrections,'Enable','On');

setappdata(0,'hMainGui',hMainGui);
fShow('Tracks');



function SaveCorrections(hMainGui)
Drift=getappdata(hMainGui.fig,'Drift'); %#ok<NASGU>
[FileName, PathName] = uiputfile({'*.mat','MAT-files (*.mat)'},'Save FIESTA Reference Transformations',fShared('GetSaveDir'));
if FileName~=0
    fShared('SetSaveDir',PathName);
    file = strtok(FileName,'.');
    save([PathName file '.mat'],'Drift');  
end

