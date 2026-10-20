function [Region,SourceFiles] = fImportImageJRegions(folder,imageSize,ScanSize)
%FIMPORTIMAGEJREGIONS Convert all foreground line ROIs below a folder.
% Returns the same Region fields as Export as Region, in pixel units.
% Conversion is completed before the caller saves or changes the GUI.
validateattributes(imageSize,{'numeric'},{'real','vector','numel',2,'integer','positive','finite'});
if ~isfolder(folder)
    error('FIESTA:InvalidROIFolder','ROI folder does not exist: %s',folder);
end
files=dir(fullfile(folder,'**','*'));
files=files(~[files.isdir]);
keep=false(size(files));
for k=1:numel(files)
    [~,name,ext]=fileparts(files(k).name);
    keep(k)=strcmpi(ext,'.roi') && ~contains(lower(name),'_bg') ...
        && ~strncmp(name,'._',2);
end
files=files(keep);
if isempty(files)
    error('FIESTA:NoROIs','No foreground .roi files were found in %s',folder);
end
SourceFiles=arrayfun(@(f) fullfile(f.folder,f.name),files,'UniformOutput',false);
% Natural order keeps ROI 2 before ROI 10 and makes numbering repeatable.
keys=cell(size(SourceFiles));
for k=1:numel(SourceFiles)
    parts=regexp(lower(SourceFiles{k}),'\d+|\D+','match');
    for j=1:numel(parts)
        if all(isstrprop(parts{j},'digit'))
            parts{j}=sprintf('%020.0f',str2double(parts{j}));
        end
    end
    keys{k}=[parts{:}];
end
[~,order]=sort(keys);
SourceFiles=SourceFiles(order);
template=struct('Area',[],'Check',1,'X',[],'Y',[],'color',[],'ScanData',[]);
Region=repmat(template,1,numel(SourceFiles));
colors=lines(24);
for k=1:numel(SourceFiles)
    [X,Y]=fReadImageJLineROI(SourceFiles{k});
    Scan=fScanGeometry(X,Y,ScanSize);
    Region(k).X=[Scan.lx fliplr(Scan.ux) Scan.lx(1)];
    Region(k).Y=[Scan.ly fliplr(Scan.uy) Scan.ly(1)];
    Region(k).Area=roipoly(imageSize(1),imageSize(2),Region(k).X,Region(k).Y);
    if ~any(Region(k).Area(:))
        error('FIESTA:ROIOutsideImage', ...
            'ROI does not overlap the loaded image. Check image dimensions/scale: %s',SourceFiles{k});
    end
    Region(k).color=colors(mod(k-1,size(colors,1))+1,:);
    Region(k).ScanData=struct('X',Scan.X,'Y',Scan.Y, ...
        'ScanSize',ScanSize,'MeasureLength',Scan.InterpD);
end
