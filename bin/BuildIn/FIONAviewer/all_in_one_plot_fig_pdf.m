function all_in_one_plot_fig_pdf(Molecule, inputArg1)
% Plot all stats from a selection of MATLAB traces as figures.
% Options add the ability to make a png or pdf plot

% BEWARE CHANNELS GET EASILY MIXED UP!! MAKE SURE TO LABEL CORRECTLY CH1 OR
% CH2 WITH YOUR MEASUREMENT BEFORE SAVING!!!

% For purposes of plotting, Ch1 is the first channel and Ch2 is the
% following neighbor
% For purposes of coding, all subscript 1's go with Ch2 and all subscript
% 2's go with Ch1. I mixed them up at some point.

offset = -15;
dx = -18; %dx = 75; %assume on-axis, horizontal MT
dy = -27; %dy = 35; %assume off-axis, vertical MT
% Note, can instead use function offset_by_law_of_means(dir, fnames) below

% Options
grid_spacing = 16;
%number of points external to mean/smooth with on both ends (so actually doubled).
% Thus, a sraw = 2 implies you are using 5 points, the original point and two on either side.
% For no running average change set to 0.
smooth_running_average_wdw = 0;


if nargin < 2
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

% option to set offsets for plots!
% [dx, dy] = offset_by_law_of_means(dir, fnames);
dx = 0;
dy = 0;
offset = -dx+40; %arbitrary value so that on and off axes don't overlap

spatiotemporal_info = zeros(2,2,1);
efo1 = [];
eco1 = [];
efo2 = [];
eco2 = [];

% loop through each selected file, assuming they are all in the same
% directory
    
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
    
    fname(ii:jj)
    fname(gg:j)

    for k = length(Molecule):-1:1
        % Have it go backwards so that it can only find the total fname.
        % For example. if you were trying to find molecule 627, but there
        % was also a molecule 66273, it would find the higher molecule
        % first. But then upon crossing to 627, it would update m
        % accordingly.

        % Contrastly, if you are trying to find molecule 66273, it will not
        % update when it reaches 627 since the str is not contained in
        % Molecule 627 name. So you are safe. Genius move!
        if ~isempty(strfind(Molecule(k).Name,fname(ii:jj)))
            % Molecule(k).Name
            % length(Molecule(k).Name)
            % strfind(Molecule(k).Name, fname(ii:jj))
            m = k;
        end
        if ~isempty(j)
        if ~isempty(strfind(Molecule(k).Name,fname(gg:j)))
            % Molecule(k).Name
            % length(Molecule(k).Name)
            % strfind(Molecule(k).Name, fname(1:j))
            n = k;
        end
        end

    end

    if ~isempty(j) %switch it
        temp = m;
        m = n;
        n = temp;

        fname2 = fname;

        strstart = strfind(fname,'nbh');
        strend = strfind(fname,'fiona');
        fname1 = strcat(fname(1:strstart-1),fname(strend:end));
        
        fprintf(strcat("Channel ", num2str(Molecule(n).Channel), ": ", Molecule(n).Name, '\n'))
    else
        fname1 = fname;
    end
    
    fprintf(strcat("Channel ", num2str(Molecule(m).Channel), ": ", Molecule(m).Name, '\n'))

    % Now dissect filename connections
    totdata = load(fullfile(dir,fname1));
    ch1data = totdata.data;
    
    efo1 = [efo1; Molecule(m).Results(:,7)];
    eco1 = [eco1; Molecule(m).Results(:,8)];

    % % Make a figure to show the efo and eco
    % figefco = figure();
    % subplot(1+~isempty(j),2,1)
    % histogram(Molecule(m).Results(:,7))
    % title(strcat(fname(1:j), " efo"))
    % legend(num2str(median(Molecule(m).Results(:,7))))
    % subplot(1+~isempty(j),2,2)
    % histogram(Molecule(m).Results(:,8))
    % title(strcat(fname(1:j), " cnt"))
    % legend(num2str(median(Molecule(m).Results(:,8))))
    % % MOLECULE

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
    
    % % magenta fit with blue data
    % % plot(time1-time1(1), ch1data.trace(idx1,1)+offset,'Color','blue','LineWidth',1.)
    % plot(time1-time1(1), ch1trace+offset,'Color','blue','LineWidth',1.)
    % plot(time1-time1(1), ch1data.trace(idx1,3)+offset,'m','LineWidth',2)
    % 
    % % plot(time1-time1(1), ch1data.trace_yx(idx1,1),'Color',[0 0 0.7],'LineWidth',1.)
    % plot(time1-time1(1), -ch1trace_yx,'Color',[0 0 0.7],'LineWidth',1.)
    % plot(time1-time1(1), -ch1data.trace_yx(idx1,3),'m','LineWidth',2)
    
    % green data slightly transparent with green line fit
    plot(time1-time1(1), ch1data.trace(idx1,1)+offset,'Color',[0 0.7 0 0.45],'LineWidth',0.4)
    plot(time1-time1(1), ch1data.trace(idx1,3)+offset,'Color',[0 0.7 0],'LineWidth',2)
    
    % Negative for now because the scipt for finding a neighbor seems to
    % have a negative cross-product result. If fixed, change this sign.
    plot(time1-time1(1), -ch1data.trace_yx(idx1,1),'Color',[0 0.4 0 0.45],'LineWidth',0.4)
    plot(time1-time1(1), -ch1data.trace_yx(idx1,3),'Color',[0 0.4 0],'LineWidth',2)

    ax = gca;
    ax.XGrid = 'off';
    ax.YGrid = 'on';
    yticks(horzcat(fliplr(-grid_spacing:-grid_spacing:ax.YLim(1)), 0:grid_spacing:ax.YLim(2)))


    if ~isempty(j) % now do the neighbor

        % try

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

            % blue fit with more transparent blue data
            plot(time2-time1(1), ch2data.trace(idx2,1)+offset+dx,'Color',[0 0 1 0.45],'LineWidth',0.4)
            plot(time2-time1(1), ch2data.trace(idx2,3)+offset+dx,'blue','LineWidth',2)

            plot(time2-time1(1), ch2data.trace_yx(idx2,1)+dy,'Color',[0 0 0.7 0.45],'LineWidth',0.4)
            plot(time2-time1(1), ch2data.trace_yx(idx2,3)+dy,'blue','LineWidth',2)
            
            % MOLECULE
            efo2 = [efo2; Molecule(n).Results(:,7)];
            eco2 = [eco2; Molecule(n).Results(:,8)];
            
            figefco = figure();
            subplot(1+~isempty(j),2,1)
            histogram(Molecule(m).Results(:,7))
            title(strcat(fname(ii:jj), " efo"))
            legend(num2str(median(Molecule(m).Results(:,7))))
            subplot(1+~isempty(j),2,2)
            histogram(Molecule(m).Results(:,8))
            title(strcat(fname(ii:jj), " cnt"))
            legend(num2str(median(Molecule(m).Results(:,8))))
            subplot(2,2,3)
            histogram(Molecule(n).Results(:,7))
            title(strcat(fname(1:j), " efo"))
            legend(num2str(median(Molecule(n).Results(:,7))))
            subplot(2,2,4)
            histogram(Molecule(n).Results(:,8))
            title(strcat(fname(1:j), " cnt"))
            legend(num2str(median(Molecule(n).Results(:,8))))
            % MOLECULE

        % catch
        %     fprintf("no associated neighbor file ... continuing with just a single channel \n")
        % 
        %     idx1 = [NaN];
        % 
        % end

    else
        fprintf("no associated neighbor file ... continuing with just a single channel \n")

        idx2 = [NaN];

        figefco = figure();
        subplot(1,2,1)
        histogram(Molecule(m).Results(:,7))
        title(strcat(fname(ii:jj), " efo"))
        legend(num2str(median(Molecule(m).Results(:,7))))
        subplot(1,2,2)
        histogram(Molecule(m).Results(:,8))
        title(strcat(fname(ii:jj), " cnt"))
        legend(num2str(median(Molecule(m).Results(:,8))))


    end
    
    st = zeros(2,2);

    if sum(~isnan(idx2)) > 0

        t1 = max([time1(1),time2(1)]);
        t2 = min([time1(end),time2(end)]);

        midx1a = find(time1 < t2); midx1b = find(time1 > t1);
        midx2a = find(time2 < t2); midx2b = find(time2 > t1);
        idx1 = idx1(intersect(midx1a,midx1b));
        idx2 = idx2(intersect(midx2a,midx2b));

        y1 = (ch1data.trace_yx(idx1(2:end),1) - ch1data.trace_yx(idx1(1:end-1),1)).^2;
        x1 = (ch1data.trace(idx1(2:end),1) - ch1data.trace(idx1(1:end-1),1) ).^2;

        st(1,1) = mean( sqrt(x1 + y1)) / sqrt(2);
        st(1,2) = mean(time1(2:end)-time1(1:end-1))*1000;

        fprintf(strcat("ch1 dx: ", num2str( round(st(1,1),2) ), " nm \n"))
        fprintf(strcat("ch1 dt: ", num2str( round(st(1,2),3) ), " ms \n"))

        y2 = (ch2data.trace_yx(idx2(2:end),1) - ch2data.trace_yx(idx2(1:end-1),1)).^2;
        x2 = (ch2data.trace(idx2(2:end),1) - ch2data.trace(idx2(1:end-1),1) ).^2;
        
        st(2,1) = mean( sqrt(x2 + y2)) / sqrt(2);
        st(2,2) = mean(time2(2:end)-time2(1:end-1))*1000;
    
        fprintf(strcat("ch2 dx: ", num2str( round(st(2,1),2)), " nm \n"))
        fprintf(strcat("ch2 deltat: ", num2str( round(st(2,2),3) ), " ms \n"))

        xycorrect = -[mean(Molecule(m).Results(:,3)) - mean(Molecule(n).Results(:,3)), mean(Molecule(m).Results(:,4)) - mean(Molecule(n).Results(:,4))];
        fprintf(strcat('Mean x, y: (', num2str(round(xycorrect(1),3)) ,', ', num2str(round(xycorrect(2),3)), ')\n'))
    
    else

        y1 = (ch1data.trace_yx(idx1(2:end),1) - ch1data.trace_yx(idx1(1:end-1),1)).^2;
        x1 = (ch1data.trace(idx1(2:end),1) - ch1data.trace(idx1(1:end-1),1) ).^2;
        
        st(1,1) = mean( sqrt(x1 + y1)) / sqrt(2);
        st(1,2) = mean(time1(2:end)-time1(1:end-1))*1000;
    
        fprintf(strcat("ch1 dx: ", num2str( round(st(1,1),2)), " nm \n"))
        fprintf(strcat("ch1 deltat: ", num2str( round(st(1,2),3) ), " ms \n"))


    end

    spatiotemporal_info(:,:,i) = st;
    
    set(0, 'currentFigure', figgs)
    fname(1:end-4)
    title(fname(1:end-4),'Interpreter','none') % why is this not working?
    xlabel("Time (s)")
    ylabel("Distance (nm)")

    saveas(figgs, strcat(dir,fname(1:end-4),"_plot"))
    saveas(figefco, strcat(dir,fname(1:end-4),"_efco"))
    
    if write_pdf > 0
        saveas(figgs, strcat(dir,fname(1:end-4),'.png'))
        fprintf("Wrote " + strcat(fname(1:end-4)) + '.png' + '\n')
        saveas(figgs, strcat(dir,fname(1:end-4),"_plot",'.pdf'))
        fprintf("Wrote " + strcat(fname(1:end-4),"_plot") + '.pdf' + '\n')   
    end
    
    close(figgs)
    close(figefco)

end

mean(spatiotemporal_info,3)

% MOLECULE
% Make a figure to show the efo and eco
figtotefco = figure();
subplot(1+~isempty(j),2,1)
histogram(efo1)
title("Ch1 efo")
legend(num2str(median(efo1)))
subplot(1+~isempty(j),2,2)
histogram(eco1)
title("Ch1 cnt")
legend(num2str(median(eco1)))
if ~isempty(efo2)
    subplot(1+~isempty(j),2,3)
    histogram(efo2)
    title("Ch2 efo")
    legend(num2str(median(efo2)))
    subplot(1+~isempty(j),2,4)
    histogram(eco2)
    title("Ch2 cnt")
    legend(num2str(median(eco2)))
end

saveas(figtotefco, strcat(dir,"efco"))
% MOLECULE


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