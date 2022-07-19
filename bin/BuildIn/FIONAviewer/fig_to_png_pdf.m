function fig_to_png_pdf(write_pdf)
%JS 220204
% MATLAB function to write all fig files selected to png with the option to
% also do pdf. Is useful for trying to present rather than having to open
% and save every individual figure

if nargin < 1
    write_pdf = 0;
else
    write_pdf = 1;
end

[fnames, dir] = uigetfile({'*.fig'},'MultiSelect','on');
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

    figgs = openfig(strcat(dir,fname));

    saveas(figgs, strcat(dir,fname(1:end-4),'.png'))
    
    if write_pdf > 0
    saveas(figgs, strcat(dir,fname(1:end-4),'.pdf'))
    fprintf("Wrote " + strcat(fname(1:end-4)) + '.pdf' + '\n')   
    end
    
    fprintf("Wrote " + strcat(fname(1:end-4)) + '.png' + '\n')   
    
    close(figgs)

end

end

