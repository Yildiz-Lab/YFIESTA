function xy_yx()
%JS 220204
% MATLAB function to switch fiona files x and y coordinates to then be
% loaded into FIONA

[fnames, dir] = uigetfile({'*.txt'},'MultiSelect','on');
%Multiselect allows user to grab more than one file and parse accordingly

% if only one file is selected, need to package into a cell array so that
% the for loop will work
if class(fnames) == 'char'
    fnames = {fnames};
end

% loop through each selected file, assuming they are all in the same
% directory
    
for i = 1:length(fnames)

    fname = fnames{i};

    xy = importdata(strcat(dir,fname));
    yx = [xy(:,2), xy(:,1)];

    dlmwrite(strcat(dir,fname(1:end-9),'yx_',fname(end-8:end)), yx, 'delimiter',' ')
    fprintf("Wrote " + strcat(fname(1:end-9),'yx_',fname(end-8:end)) + '\n')   

end

end

