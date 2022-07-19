function FileName = load_smooth_tif(filename)

%   Description:
%   Load in a tiff movie even if it is .ome which has been a problem in the past
%   IN PROGRESS: 
%       and smooth every page individually based on gaussian filtering

%   Parameters:
%       filename: (path) the full path where to find the original file
%       smoothing_sigma: (float, optional) the standard deviation by which to
%           smooth

%   Returns:
%       FileName: (string) the name of the new .tif or .tiff file now
%           formatted appropriately

%    Joseph Slivka
%    2022/05/27 
        

%[select_file, path] = uigetfile({'*.tif';'*.tiff'},'Select tif file');
%filename = fullfile(path, select_file);
% if no file selected, don't spit out an error just end function
% if select_file == 0
%     return
% end


%% Set savename
% we want to remove the .ome or .tiff option as to not confuse FIESTA so if that
% exists we will write it without it. So only save as .tif

if contains(filename, '.f') %if this is already preprocessed, don't make a new file
   [~,f,type] = fileparts(filename);
   FileName = strcat(f,type);
   fprintf("Reloading an already converted file \n If you want a different sigma, delete and reprocess \n")
   return

else
    % set to usual smoothing factor if user doesn't pass in a different one
    smoothing_sigma = str2double(fInputDlg('Define Sigma to smooth image by gaussian filtering (Default 0 is no smoothing):','0'));  
    
    k = strfind(filename, '.');

    if length(k) > 1
        k = k(1);
    end

    savename = strcat(filename(1:k-1), '.f.tif');
    [~,f,type] = fileparts(savename);
    FileName = strcat(f,type);

    if isfile(savename)
        warning('File already exists, delete file and then run')
        return
    end


    %% Load in tiff into stack
    % source for loading and saving tiffs
    %  https://www.mathworks.com/matlabcentral/answers/105739-how-to-show-tiff-stacks

    tiff_infos = imfinfo(filename); % return tiff structure, one element per image
    tiff_stack = imread(filename, 1); % read in first image

    if smoothing_sigma > 0 % do gaussian filter
    tiff_stack = imgaussfilt(tiff_stack, smoothing_sigma);
    end

    %% concatenate each successive tiff to tiff_stack
    for ii = 2 : size(tiff_infos, 1)
        temp_tiff = imread(filename, ii);
        if smoothing_sigma > 0 %do gaussian filter
        temp_tiff = imgaussfilt(temp_tiff, smoothing_sigma);% do gauss filter
        end
        tiff_stack = cat(3 , tiff_stack, temp_tiff);
    end

    %% Save tiff
    %write a Tiff file, appending each image as a new page
    for ii = 1 : size(tiff_stack, 3)
        imwrite(tiff_stack(:,:,ii) , savename , 'WriteMode' , 'append', 'Compression', 'none');
        % make sure there is no compression or else Fiesta will not load properly
    end

end

end

