function save_fit_to_mat(fullfilename)
% Joseph made a mistake fitting data but not overwriting, but he was smart
% enough to save all the figures. Therefore, he had to write a function
% that did the reverse to make up for it.


if nargin < 1
    % [filename,pathname] = uigetfile('MultiSelect','on')
    [filename,pathname] = uigetfile()
    fullfilename = fullfile(pathname,filename);
end

[d,f,e] = fileparts(fullfilename);
% if contains(f,'._') %JS Edit to delete extra ._ from an error
% f = f(3:end);
% end

if contains(e,'fig')
    % load everything from the figure
    fig = openfig(fullfilename);
    ax = get(fig,'Children');
    xfitline = findobj(ax,'DisplayName','data1'); xfit = xfitline.YData;
    yfitline = findobj(ax,'DisplayName','data2'); yfit = yfitline.YData;
    
    savematfile = fullfile(d,strcat(f,'.mat'));
    steptrace = load(savematfile);
    
    data = steptrace.data;
    data.trace(:,3) = xfit;
    data.trace(:,4) = yfit;
    data.trace(:,5) = zeros(1,length(xfit));
    steps = xfit(2:end)-xfit(1:end-1);
    data.trace(2:end,5) = abs(steps) > 0;
    
    
    % % also populate trace_yx
    % yxfullfilename = fullfile(d,strcat(f,'_yx',e))
    % fig2 = openfig(yxfullfilename);
    % ax2 = get(fig2,'Children');
    % xfitline = findobj(ax2,'DisplayName','data1'); xfit = xfitline.YData;
    % yfitline = findobj(ax2,'DisplayName','data2'); yfit = yfitline.YData;
    % 
    % data.trace_yx(:,1) = data.trace(:,2);
    % data.trace_yx(:,2) = data.trace(:,1);
    % data.trace_yx(:,3) = xfit;
    % data.trace_yx(:,4) = yfit;
    % data.trace_yx(:,5) = zeros(1,length(xfit));
    % data.trace_yx(:,6) = zeros(1,length(xfit));
    % steps = xfit(2:end)-xfit(1:end-1);
    % data.trace_yx(2:end,5) = abs(steps) > 0;

    save(savematfile,'data');

end


end