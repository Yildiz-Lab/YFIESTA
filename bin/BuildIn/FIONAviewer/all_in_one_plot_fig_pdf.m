function all_in_one_plot_fig_pdf(inputArg1)
% Plot all stats from a selection of MATLAB traces as figures.
% Options add the ability to make a png or pdf plot

offset = 35;
dx = 0; %assume on-axis, horizontal MT
dy = 0; %assume off-axis, vertical MT
grid_spacing = 16;

if nargin < 1
    write_pdf = 0;
else
    write_pdf = 1;
end

[fnames, dir] = uigetfile({'*.mat'},'MultiSelect','on');
%Multiselect allows user to grab more than one file and parse accordingly

% if only one file is selected, need to package into a cell array so that
% the for loop will work
if class(fnames) == 'char'
    fnames = {fnames};
end

% loop through each selected file, assuming they are all in the same
% directory

spatiotemporal_info = zeros(2,2,1);
    
for i = 1:length(fnames)
    try
        fname = fnames{i};
    catch
        fprintf("Finished combining xy and yx traces. Aborting... \n")
        return
    end
    
    % Now dissect filename connections

    totdata = load(fullfile(dir,fname));
    ch2data = totdata.data;

    if isfield(ch2data, 'time')
        time2 = ch2data.time;
        time2 = ch2data.time(~isnan(ch2data.time));
    else
        time2 = 1:size(ch2data.xy,1);
    end
    
    figgs = figure();
    hold on

    idx2 = find(~isnan(ch2data.trace(:,1)));
    
    % plot(time2-time2(1), ch2data.trace(idx2,1)+offset,'Color',[0 0.7 0],'LineWidth',1.)
    % plot(time2-time2(1), ch2data.trace(idx2,3)+offset,'k','LineWidth',2)
    % 
    % plot(time2-time2(1), ch2data.trace_yx(idx2,1),'Color',[0 0.4 0],'LineWidth',1.)
    % plot(time2-time2(1), ch2data.trace_yx(idx2,3),'k','LineWidth',2)

    plot(time2-time2(1), ch2data.trace(idx2,1)+offset,'Color','blue','LineWidth',1.)
    plot(time2-time2(1), ch2data.trace(idx2,3)+offset,'m','LineWidth',2)

    plot(time2-time2(1), ch2data.trace_yx(idx2,1),'Color',[0 0 0.7],'LineWidth',1.)
    plot(time2-time2(1), ch2data.trace_yx(idx2,3),'m','LineWidth',2)

    ax = gca;
    ax.XGrid = 'off';
    ax.YGrid = 'on';
    yticks(horzcat(fliplr(-grid_spacing:-grid_spacing:ax.YLim(1)), 0:grid_spacing:ax.YLim(2)))
    
    y2 = (ch2data.trace_yx(idx2(2:end),1) - ch2data.trace_yx(idx2(1:end-1),1)).^2;
    x2 = (ch2data.trace(idx2(2:end),1) - ch2data.trace(idx2(1:end-1),1) ).^2;
    
    st = zeros(2,2);
    st(2,1) = mean( sqrt(x2 + y2)) / sqrt(2);
    st(2,2) = mean(time2(2:end)-time2(1:end-1))*1000;

    fprintf(strcat("ch2 dx: ", num2str( round(st(2,1),2)), " nm \n"))
    fprintf(strcat("ch2 deltat: ", num2str( round(st(2,2),3) ), " ms \n"))


    if contains(fname,'nbh')
        strstart = strfind(fname,'nbh');
        strend = strfind(fname,'fiona');
        fname1 = strcat(fname(1:strstart-1),fname(strend:end));
        % fprintf("nbh")
        
        % try
            totdata = load(fullfile(dir,fname1));
            ch1data = totdata.data;
    
            if isfield(ch1data, 'time')
                time1 = ch1data.time;
                time1 = ch1data.time(~isnan(ch1data.time));
            else
                time1 = 1:size(ch1data.xy,1);
            end
    
            idx1 = find(~isnan(ch1data.trace(:,1)));
            % time2 is not a typo, want to make sure you are normalizing to the
            % same time
            plot(time1-time2(1), ch1data.trace(idx1,1)+offset+dx,'Color',[0 0.7 0],'LineWidth',1.)
            plot(time1-time2(1), ch1data.trace(idx1,3)+offset+dy,'k','LineWidth',2)
    
            plot(time1-time2(1), ch1data.trace_yx(idx1,1),'Color',[0 0.4 0],'LineWidth',1.)
            plot(time1-time2(1), ch1data.trace_yx(idx1,3),'k','LineWidth',2)
            
            y2 = (ch1data.trace_yx(idx1(2:end),1) - ch1data.trace_yx(idx1(1:end-1),1)).^2;
            x2 = (ch1data.trace(idx1(2:end),1) - ch1data.trace(idx1(1:end-1),1) ).^2;

            st(1,1) = mean( sqrt(x2 + y2)) / sqrt(2);
            st(1,2) = mean(time1(2:end)-time1(1:end-1))*1000;

            fprintf(strcat("ch1 dx: ", num2str( round(st(1,1),2) ), " nm \n"))
            fprintf(strcat("ch1 dt: ", num2str( round(st(1,2),3) ), " ms \n"))

        % catch
            % fprintf("no associated neighbor file ... continuing with just a single channel \n")

        % end


    end

    title(fname(1:end-4),'Interpreter','none')
    xlabel("Time (s)")
    ylabel("Distance (nm)")
    
    spatiotemporal_info(:,:,i) = st;

    saveas(figgs, strcat(dir,fname(1:end-4),"_plot"))

    % saveas(figgs, strcat(dir,fname(1:end-4),'.png'))
    % fprintf("Wrote " + strcat(fname(1:end-4)) + '.png' + '\n') \
    
    if write_pdf > 0
        saveas(figgs, strcat(dir,fname(1:end-4),'.png'))
        fprintf("Wrote " + strcat(fname(1:end-4)) + '.png' + '\n')
        saveas(figgs, strcat(dir,fname(1:end-4),"_plot",'.pdf'))
        fprintf("Wrote " + strcat(fname(1:end-4),"_plot") + '.pdf' + '\n')   
    end
    
    close(figgs)

end

mean(spatiotemporal_info,3)