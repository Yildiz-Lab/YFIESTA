function all_in_one_plot_fig_pdf(Molecule, inputArg1)
% Plot all stats from a selection of MATLAB traces as figures.
% Options add the ability to make a png or pdf plot

% BEWARE CHANNELS GET EASILY MIXED UP!! MAKE SURE TO LABEL CORRECTLY CH1 OR
% CH2 WITH YOUR MEASUREMENT BEFORE SAVING!!!

offset = -15;
dx = -18; %dx = 75; %assume on-axis, horizontal MT
dy = -27; %dy = 35; %assume off-axis, vertical MT
% Note, can instead use function offset_by_law_of_means(dir, fnames) below

grid_spacing = 16;

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
[dx, dy] = offset_by_law_of_means(dir, fnames);
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
    
    % Now dissect filename connections

    totdata = load(fullfile(dir,fname));
    ch2data = totdata.data;

    % MOLECULE
    % find the molecule in the data Molecule if loaded
    jj = strfind(fname,'_fiona')-1; %nbh or initial
    j = strfind(fname,'_nbh')-1; %if there is a neighbor
    ii = 1;
    if ~isempty(j)
        ii = j+6;
    end
    
    fname(ii:jj)
    fname(1:j)

    for k = 1:length(Molecule)
        if ~isempty(strfind(Molecule(k).Name,fname(ii:jj)))
            m = k;
        end
        if ~isempty(j)
        if ~isempty(strfind(Molecule(k).Name,fname(1:j)))
            n = k;
        end
        end
    end

    if ~isempty(j) %switch it
        temp = m;
        m = n;
        n = temp;
    end
    
    fprintf('Channel 2 \n')
    Molecule(m).Name
    Molecule(m).Channel
    efo2 = [efo2; Molecule(m).Results(:,7)];
    eco2 = [eco2; Molecule(m).Results(:,8)];

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
            plot(time1-time2(1), ch1data.trace(idx1,3)+offset+dx,'k','LineWidth',2)
    
            plot(time1-time2(1), ch1data.trace_yx(idx1,1)+dy,'Color',[0 0.4 0],'LineWidth',1.)
            plot(time1-time2(1), ch1data.trace_yx(idx1,3)+dy,'k','LineWidth',2)
            
            % MOLECULE
            efo1 = [efo1; Molecule(n).Results(:,7)];
            eco1 = [eco1; Molecule(n).Results(:,8)];
            
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

        idx1 = [NaN];


    end
    
    st = zeros(2,2);

    if sum(~isnan(idx1)) > 0

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
    
    else

        y2 = (ch2data.trace_yx(idx2(2:end),1) - ch2data.trace_yx(idx2(1:end-1),1)).^2;
        x2 = (ch2data.trace(idx2(2:end),1) - ch2data.trace(idx2(1:end-1),1) ).^2;
        
        st(2,1) = mean( sqrt(x2 + y2)) / sqrt(2);
        st(2,2) = mean(time2(2:end)-time2(1:end-1))*1000;
    
        fprintf(strcat("ch2 dx: ", num2str( round(st(2,1),2)), " nm \n"))
        fprintf(strcat("ch2 deltat: ", num2str( round(st(2,2),3) ), " ms \n"))


    end

    spatiotemporal_info(:,:,i) = st;
    
    set(0, 'currentFigure', figgs)
    fname(1:end-4)
    title(fname(1:end-4),'Interpreter','none') % why is this not working?
    xlabel("Time (s)")
    ylabel("Distance (nm)")

    saveas(figgs, strcat(dir,fname(1:end-4),"_plot"))
    saveas(figefco, strcat(dir,fname(1:end-4),"_efco"))

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

% MOLECULE
% Make a figure to show the efo and eco
figtotefco = figure();
subplot(1+~isempty(j),2,1)
histogram(efo2)
title("Ch2 efo")
legend(num2str(median(efo2)))
subplot(1+~isempty(j),2,2)
histogram(eco2)
title("Ch2 cnt")
legend(num2str(median(eco2)))
if ~isempty(efo1)
    subplot(2,2,3)
    histogram(efo1)
    title("Ch1 efo")
    legend(num2str(median(efo1)))
    subplot(2,2,4)
    histogram(eco1)
    title("Ch1 cnt")
    legend(num2str(median(eco1)))
end

saveas(figtotefco, strcat(dir,"efco"))
% MOLECULE