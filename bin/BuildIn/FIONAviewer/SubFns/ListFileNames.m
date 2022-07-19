function fileList = ListFileNames( dirName )
%LISTFILES list all .bin files in directory

D = dir([dirName,'\*.bin']);

fileList = {};

for i=1:length(D)
    fileList = [fileList; D(i).name]; %#ok<AGROW>
end

fileList = char(fileList);