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
    
    % MOLECULE
    % find the molecule in the data Molecule if loaded
    jj = strfind(fname,'_fiona')-1; %nbh or initial
    j = strfind(fname,'_nbh')-1; %if there is a neighbor
    ii = 1; gg = 1;
    if ~isempty(j)
        ii = j+6;
    end

    if contains(fname(1:j),'_') % actually a more complicated name
        gg = strfind(fname(1:j),'_')+1;
        if length(gg) > 1
        gg = gg(end);
        end
    end
    
    % fname(ii:jj)
    % fname(gg:j)

    if ~isempty(j) %switch it

        fname2 = fname;

        strstart = strfind(fname,'nbh');
        strend = strfind(fname,'fiona');
        % JS Edit works for filename types of
        % (mol#)_nbh_(mol#)_(whatever_whatever)_fiona.mat
        aa = strfind(fname,'_');
        aa = aa(aa > strstart);
        strend = aa(2)+1; %Pull second since first is to prelude the nbh number, then following the neighbor number
        fname1 = strcat(fname(1:strstart-1),fname(strend:end));
    else
        fname1 = fname;
    end
    % fname1
    % fname2

        % Now dissect filename connections
    totdata = load(fullfile(dir,fname1));
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

    if ~isempty(j) % now do the neighbor

        totdata = load(fullfile(dir,fname2));
        ch2data = totdata.data;

        if isfield(ch2data, 'time')
            time2 = ch2data.time;
            time2 = ch2data.time(~isnan(ch2data.time));
        else
            time2 = 1:size(ch2data.xy,1);
        end

        idx2 = find(~isnan(ch2data.trace(:,1)));
        % time1 is not a typo, want to make sure you are normalizing to the
        % same time

        % Smooth by mean window
        ch2trace = smooth_running_average(ch2data.trace(idx2,1), smooth_running_average_wdw);
        ch2trace_yx = smooth_running_average(ch2data.trace_yx(idx2,1), smooth_running_average_wdw);

        % % green data with black line fit
        % % plot(time2-time1(1), ch2data.trace(idx2,1)+offset+dx,'Color',[0 0.7 0],'LineWidth',1.)
        % plot(time2-time1(1), ch2trace+offset+dx,'Color',[0 0.7 0],'LineWidth',1.)
        % plot(time2-time1(1), ch2data.trace(idx2,3)+offset+dx,'k','LineWidth',2)
        % 
        % % plot(time2-time1(1), ch2data.trace_yx(idx2,1)+dy,'Color',[0 0.4 0],'LineWidth',1.)
        % plot(time2-time1(1), ch2trace_yx+dy,'Color',[0 0.4 0],'LineWidth',1.)
        % plot(time2-time1(1), ch2data.trace_yx(idx2,3)+dy,'k','LineWidth',2)

        % % blue fit with more transparent blue data
        % plot(time2-time1(1), ch2data.trace(idx2,1)+offset,'Color',[0 0 1 0.45],'LineWidth',0.4)
        % plot(time2-time1(1), ch2data.trace(idx2,3)+offset,'blue','LineWidth',2)
        % 
        % plot(time2-time1(1), ch2data.trace_yx(idx2,1),'Color',[0 0 0.7 0.45],'LineWidth',0.4)
        % plot(time2-time1(1), ch2data.trace_yx(idx2,3),'blue','LineWidth',2)

        % blue data with yellow line fit
        plot(time2-time1(1), ch2trace+offset,'Color','b','LineWidth',1.)
        plot(time2-time1(1), ch2data.trace(idx2,3)+offset,'Color',[212,175,55]/255,'LineWidth',2)

        plot(time2-time1(1), ch2trace_yx,'Color','b','LineWidth',1.)
        plot(time2-time1(1), ch2data.trace_yx(idx2,3),'Color',[212,175,55]/255,'LineWidth',2)
    end


    ax = gca;
    ax.XGrid = 'off';
    ax.YGrid = 'on';
    yticks(horzcat(fliplr(-grid_spacing:-grid_spacing:ax.YLim(1)), 0:grid_spacing:ax.YLim(2)))

    title(fname(1:end-4),'Interpreter','none') % why is this not working?
    xlabel("Time (s)")
    ylabel("Distance (nm)")

    set(ax, ...
        'FontName', 'Arial', ...
        'FontSize', 10, ...
        'TickDir', 'out', ...
        'LineWidth', 1, ...
        'Box', 'off', ...
        'XColor', 'k', ...
        'YColor', 'k');
    
    % Optional: Update label font/colors if needed
    if ~isempty(get(ax, 'XLabel'))
        set(get(ax, 'XLabel'), 'FontName', 'Arial', 'FontSize', 10, 'Color', 'k');
    end
    if ~isempty(get(ax, 'YLabel'))
        set(get(ax, 'YLabel'), 'FontName', 'Arial', 'FontSize', 10, 'Color', 'k');
    end

% Optional: white background
set(gcf, 'Color', 'w');

    saveas(figgs, strcat(dir,fname(1:end-4),"_plot"))
    
    if write_pdf > 0
        saveas(figgs, strcat(dir,fname(1:end-4),'.png'))
        fprintf("Wrote " + strcat(fname(1:end-4)) + '.png' + '\n')
        saveas(figgs, strcat(dir,fname(1:end-4),"_plot",'.pdf'))
        fprintf("Wrote " + strcat(fname(1:end-4),"_plot") + '.pdf' + '\n')   
    end
    close(figgs)

    % adding 2-Color separation histograms, looking for differences in sigma
    % offset_by_raw_means_v2(fullfile(dir,fname1))
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