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

refChannel = str2double(fInputDlg('Enter reference channel number (stays fixed):','1'));
if isnan(refChannel)
    return;
end
adjChannel = str2double(fInputDlg('Enter channel number to adjust:','2'));
if isnan(adjChannel)
    return;
end


% DC260910
% Revert any correction currently baked into this channel's tracks.
% Per-molecule, so it's safe even if only some of them are corrected.
for i = 1:length(Molecule)
    if Molecule(i).Channel == adjChannel && Molecule(i).MultiChannelCorrected
        Molecule(i).Results(:,3) = Molecule(i).Results(:,3) + Molecule(i).MultiChannelDelta(:,1);
        Molecule(i).Results(:,4) = Molecule(i).Results(:,4) + Molecule(i).MultiChannelDelta(:,2);
        Molecule(i).Results(:,6) = sqrt( (Molecule(i).Results(:,3) - Molecule(i).Results(1,3)).^2 + (Molecule(i).Results(:,4) - Molecule(i).Results(1,4)).^2 );
        Molecule(i).MultiChannelCorrected = false;
        Molecule(i).MultiChannelDelta = [];
    end
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

for i = 1:length(Molecule)
    if Molecule(i).Channel == adjChannel
        rawXY = double(Molecule(i).Results(:,3:4));
        dx = Fx(rawXY);
        dy = Fy(rawXY);
        Molecule(i).Results(:,3) = Molecule(i).Results(:,3) - dx;
        Molecule(i).Results(:,4) = Molecule(i).Results(:,4) - dy;
        Molecule(i).Results(:,6) = sqrt( (Molecule(i).Results(:,3) - Molecule(i).Results(1,3)).^2 + (Molecule(i).Results(:,4) - Molecule(i).Results(1,4)).^2 );
        Molecule(i).MultiChannelDelta = [dx, dy];
        Molecule(i).MultiChannelCorrected = true;
    end
end

setappdata(hMainGui.fig,'MaskCorrect',F);
setappdata(hMainGui.fig,'MaskCorrectFiles',{ch1file, ch2file});
setappdata(hMainGui.fig,'MaskCorrectChannels',[refChannel, adjChannel]);

set(hMainGui.Menu.mApplyMultiCorrections,'Enable','On');
set(hMainGui.Menu.mSaveMaskCorrections,'Enable','On');

setappdata(0,'hMainGui',hMainGui);
fShared('UpdateMenu',hMainGui);
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

F = getappdata(hMainGui.fig,'MaskCorrect');
channels = getappdata(hMainGui.fig,'MaskCorrectChannels');
if isempty(F) || isempty(channels)
    return;
end
Fx = F{1}; Fy = F{2};
adjChannel = channels(2);

idx = find([Molecule.Channel] == adjChannel);
if isempty(idx)
    return;
end
corrected = false(1,numel(idx));
for k = 1:numel(idx)
    val = Molecule(idx(k)).MultiChannelCorrected;
    delta = Molecule(idx(k)).MultiChannelDelta;
    % A molecule only counts as corrected if it has BOTH a true flag
    % AND a real delta to invert -- anything else is treated as
    % not-corrected rather than trusted blindly.
    corrected(k) = ~isempty(val) && val && ~isempty(delta) && size(delta,2)==2;
end

if all(corrected)
    % fully corrected -> invert all, no prompt
    for i = idx
        Molecule(i).Results(:,3) = Molecule(i).Results(:,3) + Molecule(i).MultiChannelDelta(:,1);
        Molecule(i).Results(:,4) = Molecule(i).Results(:,4) + Molecule(i).MultiChannelDelta(:,2);
        Molecule(i).Results(:,6) = sqrt( (Molecule(i).Results(:,3) - Molecule(i).Results(1,3)).^2 + (Molecule(i).Results(:,4) - Molecule(i).Results(1,4)).^2 );
        Molecule(i).MultiChannelCorrected = false;
        Molecule(i).MultiChannelDelta = [];
    end
else
    % off or mixed -> bring everyone up to corrected, touching only
    % the ones that aren't already
    for i = idx
        val = Molecule(i).MultiChannelCorrected;
        if isempty(val) || ~val
            rawXY = double(Molecule(i).Results(:,3:4));
            dx = Fx(rawXY);
            dy = Fy(rawXY);
            Molecule(i).Results(:,3) = Molecule(i).Results(:,3) - dx;
            Molecule(i).Results(:,4) = Molecule(i).Results(:,4) - dy;
            Molecule(i).Results(:,6) = sqrt( (Molecule(i).Results(:,3) - Molecule(i).Results(1,3)).^2 + (Molecule(i).Results(:,4) - Molecule(i).Results(1,4)).^2 );
            Molecule(i).MultiChannelDelta = [dx, dy];
            Molecule(i).MultiChannelCorrected = true;
        end
    end
end

fShared('UpdateMenu',hMainGui);
fShow('Tracks');
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
Fx = F{1}; Fy = F{2};

% Reconstruct true raw positions (undoing any correction currently on
% this channel, per molecule) and compute+validate the new correction
% against those raw positions -- before touching any Molecule data, so
% a bad correction never leaves tracks half-reverted.
rawPos = cell(1,length(Molecule));
deltas = cell(1,length(Molecule));
badIdx = [];
for i = 1:length(Molecule)
    if Molecule(i).Channel == adjChannel
        if Molecule(i).MultiChannelCorrected
            raw = double(Molecule(i).Results(:,3:4)) + double(Molecule(i).MultiChannelDelta);
        else
            raw = double(Molecule(i).Results(:,3:4));
        end
        rawPos{i} = raw;
        dx = Fx(raw);
        dy = Fy(raw);
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
        Molecule(i).Results(:,3) = rawPos{i}(:,1) - deltas{i}(:,1);
        Molecule(i).Results(:,4) = rawPos{i}(:,2) - deltas{i}(:,2);
        Molecule(i).Results(:,6) = sqrt( (Molecule(i).Results(:,3) - Molecule(i).Results(1,3)).^2 + (Molecule(i).Results(:,4) - Molecule(i).Results(1,4)).^2 );
        Molecule(i).MultiChannelDelta = deltas{i};
        Molecule(i).MultiChannelCorrected = true;
    end
end

setappdata(hMainGui.fig,'MaskCorrect',F);
setappdata(hMainGui.fig,'MaskCorrectChannels',channels);

set(hMainGui.Menu.mApplyMultiCorrections,'Enable','On');
set(hMainGui.Menu.mSaveMaskCorrections,'Enable','On');

setappdata(0,'hMainGui',hMainGui);
fShared('UpdateMenu',hMainGui);
fShow('Tracks');



function SaveCorrections(hMainGui)
Drift=getappdata(hMainGui.fig,'Drift'); %#ok<NASGU>
[FileName, PathName] = uiputfile({'*.mat','MAT-files (*.mat)'},'Save FIESTA Reference Transformations',fShared('GetSaveDir'));
if FileName~=0
    fShared('SetSaveDir',PathName);
    file = strtok(FileName,'.');
    save([PathName file '.mat'],'Drift');  
end

