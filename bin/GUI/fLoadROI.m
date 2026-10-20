function [Region,PathName,Alignment] = fLoadROI(hMainGui,source,alignRegions)
%FLOADROI Select/load saved regions or import an ImageJ ROI folder.
% All ImageJ-specific UI and file handling lives outside fLeftPanel/LoadRegion.
% Optional source and alignment arguments permit scripted imports. Alignment
% defaults to false on every call and never overwrites an existing MAT file.
Region=[];
PathName='';
Alignment=[];
if nargin<3, alignRegions=false; end
interactive=nargin<2;
if interactive
    [choice,alignRegions]=fRegionLoadDialog;
    if strcmp(choice,'MAT file')
        [FileName,PathName]=uigetfile({'*.mat','FIESTA Regions (*.mat)'}, ...
            'Load FIESTA Regions',fShared('GetLoadDir'));
        if isequal(FileName,0), return; end
        source=fullfile(PathName,FileName);
    elseif strcmp(choice,'ImageJ ROI folder')
        source=uigetdir(fShared('GetLoadDir'),'Select ImageJ ROI folder (all subfolders included)');
        if isequal(source,0), return; end
    else
        return;
    end
end
try
    if isfolder(source)
        [Region,PathName,Alignment]=LoadImageJRegions(hMainGui,source,alignRegions);
    else
        PathName=fileparts(source);
        Region=fLoad(source,'Region');
        if alignRegions
            [Region,Alignment]=AlignToCurrentMovie(Region,hMainGui);
        end
    end
    if ~isempty(PathName) && PathName(end)~=filesep
        PathName=[PathName filesep];
    end
catch exception
    if ~interactive, rethrow(exception); end
    Region=[];
    fMsgDlg(exception.message,'error');
end

function [Region,PathName,Alignment] = LoadImageJRegions(hMainGui,folder,alignRegions)
% Convert the entire folder before saving or adding anything to the GUI.
global Stack;
if isempty(Stack) || isempty(Stack{1})
    error('FIESTA:NoImage','Load an image before importing ImageJ ROIs.');
end
imageSize=[size(Stack{1},1) size(Stack{1},2)];
Region=fImportImageJRegions(folder,imageSize,hMainGui.Values.ScanSize);
Alignment=[];
if alignRegions
    [Region,Alignment]=AlignToCurrentMovie(Region,hMainGui);
end
% Save beside the selected folder; keep previous conversions intact.
folder=char(folder);
while numel(folder)>1 && (folder(end)=='/' || folder(end)==char(92))
    folder(end)=[];
end
[PathName,name]=fileparts(folder);
regionFile=fullfile(PathName,[name ' regions.mat']);
suffix=2;
while exist(regionFile,'file')
    regionFile=fullfile(PathName,sprintf('%s regions_%d.mat',name,suffix));
    suffix=suffix+1;
end
save(regionFile,'Region');
% Load the saved file through the original MAT reader.
Region=fLoad(regionFile,'Region');
fprintf('Imported %d ImageJ ROIs (pixels, ScanSize %g). Saved: %s\n', ...
    numel(Region),hMainGui.Values.ScanSize,regionFile);

function [Region,Alignment] = AlignToCurrentMovie(Region,hMainGui)
global Stack;
channel=1;
if isfield(hMainGui.Values,'FrameIdx')
    channel=hMainGui.Values.FrameIdx(1);
end
if isempty(Stack) || channel<1 || channel>numel(Stack) || isempty(Stack{channel})
    error('FIESTA:NoImage','Load a movie before aligning regions.');
end
Projection=max(Stack{channel},[],3);
[Region,Alignment]=fAlignRegionsToProjection(Region,Projection);
Alignment.Channel=channel;
fprintf('Aligned regions to channel %d maximum z-projection: dX=%+.3f, dY=%+.3f pixels (score %.3f -> %.3f).\n', ...
    channel,Alignment.ShiftXY,Alignment.ScoreBefore,Alignment.ScoreAfter);
