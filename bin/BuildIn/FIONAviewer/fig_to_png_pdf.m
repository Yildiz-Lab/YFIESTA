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
    try
        fname = fnames{i};
    catch
        fprintf("Finished combining xy and yx traces. Aborting... \n")
        return
    end
    % new option to grab the data from yx and put it on the original
    % figure, but only if it is included in the set selected already
    fyxname = strcat(fnames{i}(1:end-9), 'yx', fnames{i}(end-9:end));
    idx = find(contains(fnames,fyxname));
    if ~isempty(idx)
        figgyx = openfig(strcat(dir,fyxname));
        axchild = get(gca,'Children');
        yx = findobj(axchild, 'DisplayName', 'data1');
        X = yx.XData; Y = yx.YData;
        close(figgyx);
        fnames(idx) = [];
    end

    figgs = openfig(strcat(dir,fname));
    % only if yx is included should you copy it over
    if ~isempty(idx)
        axchild = get(gca,'Children');
        delete(findobj(axchild, 'DisplayName', 'data2'));
        hold on
        plot(X,Y,'Color',[0.3,0.3,0.9],'LineStyle','-','LineWidth',2,'MarkerSize',6,'Marker','none','MarkerFaceColor','None','DisplayName','data2')
    end
    
    saveas(figgs, strcat(dir,fname(1:end-4),'.png'))
    
    if write_pdf > 0
    saveas(figgs, strcat(dir,fname(1:end-4),'.pdf'))
    fprintf("Wrote " + strcat(fname(1:end-4)) + '.pdf' + '\n')   
    end
    
    fprintf("Wrote " + strcat(fname(1:end-4)) + '.png' + '\n')   
    
    close(figgs)

end

end

