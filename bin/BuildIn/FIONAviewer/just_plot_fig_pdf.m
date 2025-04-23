% plot_fig_pdf
% for 1C simple until I get around to doing it in all_in_one_plot_fig_pdf

offset = 40; %arbitrary value so that on and off axes don't overlap
grid_spacing = 16;
write_pdf = 0;
smooth_running_average_wdw = 0;

[fnames, dir] = uigetfile({'*.mat'},'MultiSelect','on');
%Multiselect allows user to grab more than one file and parse accordingly

% if only one file is selected, need to package into a cell array so that
% the for loop will work
if class(fnames) == 'char'
    fnames = {fnames};
end



for i = 1:length(fnames)
    try
        fname = fnames{i};
    catch
        fprintf("Finished combining xy and yx traces. Aborting... \n")
        return
    end

        % Now dissect filename connections
    totdata = load(fullfile(dir,fname));
    ch1data = totdata.data;

    if isfield(ch1data, 'time')
        time1 = ch1data.time;
        time1 = ch1data.time(~isnan(ch1data.time));
    else
        time1 = 1:size(ch1data.xy,1);
    end
    
    figgs = figure();
    hold on

    idx1 = find(~isnan(ch1data.trace(:,1)));

    % Smooth by mean window
    ch1trace = smooth_running_average(ch1data.trace(idx1,1), smooth_running_average_wdw);
    ch1trace_yx = smooth_running_average(ch1data.trace_yx(idx1,1), smooth_running_average_wdw);
    
    % % % green data with black line fit
    % plot(time1-time1(1), ch1data.trace(idx1,1)+offset,'Color',[0 0.7 0],'LineWidth',1.)
    % plot(time1-time1(1), ch1data.trace(idx1,3)+offset,'Color','k','LineWidth',2)
    % 
    % % Negative for now because the scipt for finding a neighbor seems to
    % % have a negative cross-product result. If fixed, change this sign.
    % plot(time1-time1(1), -ch1data.trace_yx(idx1,1),'Color',[0 0.4 0],'LineWidth',1.)
    % plot(time1-time1(1), -ch1data.trace_yx(idx1,3),'Color','k','LineWidth',2)

    % % black data with red line fit
    plot(time1-time1(1), ch1data.trace(idx1,1)+offset,'Color','k','LineWidth',1.)
    plot(time1-time1(1), ch1data.trace(idx1,3)+offset,'Color','r','LineWidth',2)

    % Negative for now because the scipt for finding a neighbor seems to
    % have a negative cross-product result. If fixed, change this sign.
    plot(time1-time1(1), -ch1data.trace_yx(idx1,1),'Color',[0.2 0.2 0.2],'LineWidth',1.)
    plot(time1-time1(1), -ch1data.trace_yx(idx1,3),'Color','r','LineWidth',2)

    ax = gca;
    ax.XGrid = 'off';
    ax.YGrid = 'on';
    yticks(horzcat(fliplr(-grid_spacing:-grid_spacing:ax.YLim(1)), 0:grid_spacing:ax.YLim(2)))

    title(fname(1:end-4),'Interpreter','none') % why is this not working?
    xlabel("Time (s)")
    ylabel("Distance (nm)")

    saveas(figgs, strcat(dir,fname(1:end-4),"_plot"))
    
    if write_pdf > 0
        saveas(figgs, strcat(dir,fname(1:end-4),'.png'))
        fprintf("Wrote " + strcat(fname(1:end-4)) + '.png' + '\n')
        saveas(figgs, strcat(dir,fname(1:end-4),"_plot",'.pdf'))
        fprintf("Wrote " + strcat(fname(1:end-4),"_plot") + '.pdf' + '\n')   
    end
    close(figgs)

end


function A = smooth_running_average(A, smooth_running_average_wdw)

if smooth_running_average_wdw < 1
    return
else

ll = length(A);
A = [nan(smooth_running_average_wdw,1); A; nan(smooth_running_average_wdw,1)];
for c = 1:ll
    A(c) = mean(A(c:c+2*smooth_running_average_wdw),'omitnan');
end
% now delete the excess
A(end+1-2*smooth_running_average_wdw:end) = [];

end
end