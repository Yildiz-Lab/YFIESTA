function [dx, dy] = offset_by_raw_means_v2(molecule_filepath)
% Use law of means to just guess the offset for a set of traces. Should
% converge by the CLT, so make sure you have a sufficient number of traces.

% Parameters
% Molecule : path to Molecule data structure(s)

% Returns
% don : (float) average difference of ch2 - ch1 in on-axis, thus correction
% to ch1
% doff : (float) average difference of ch2 - ch1 in off-axis, thus
% correction to ch1

%% We begin with the normal script for opening the data for all files and their neighbors.
% Taken verbadum from all_in_one_plot_fig_pdf

% % debugging
% [dir, fname, ~] = fileparts('F:\MINFLUX JS\Kinesin 2C\241216\241216_03\241216-110351_minflux_655_LP7p5_bg40k_555_LP4_bg15k_fiesta_annotated');
% fnames = {fname};

decimation_factor = 10;
% cutoff for std dev max or pure nm, will take the minimum
% set other to something high if you want to use 
std_max_factor = 1.0;
std_max = 25;
% define the beginning index to use for alignment, a useful trick if you have lots of tails you don't want applied to your data, especially with gold beads
bidx = 1;

if nargin < 1
    [fnames, dir] = uigetfile({'*.mat'},'Select files to calculate linear correction','MultiSelect','on')
else
    [dir, fnames, ~] = fileparts(molecule_filepath)
end
% if only one file is selected, need to package into a cell array so that
% the for loop will work
if class(fnames) == 'char'
    fnames = {fnames};
end

% load the molecule_datafile cell

dx_all = [];
dy_all = [];

molinfo = [];
molnum = [];

for i = 1:length(fnames)
    try
        fname = fnames{i};
    catch
        fprintf("Finished combining xy and yx traces. Aborting... \n")
        return
    end
    
    % Now dissect filename connections

    data = load(fullfile(dir,fname));
    % if isfield(data,'Molecule')
        Molecule = data.Molecule;
    % else
    %     Molecule(1).Channel = 1;
    %     Molecule(1).Results = nan(length(data.trace(:,1)),4);
    % 
    %     Molecule(1).Results(:,2) = 
    % end

    for j = 1:length(Molecule)
        % check it is in Channel 1
        if Molecule(j).Channel == 1
            % Now look for temporally overlapping neighbors. We can
            % minimize our search to just the next couple molecules as this
            % would be the only places the molecules would exist
            [ch1tmin, ch1tmax] = bounds(Molecule(j).Results(:,2));
            ch1_time_footprint = [ch1tmin, ch1tmax];
            
            checknext = j+min(3,length(Molecule)-j); %only need to check the next 3
            
            for k = j+1:checknext
                
                % check that there is a temporal overlap before we start
                % matching the molecules as temporal neighbors
                if ~( all( ch1_time_footprint < Molecule(k).Results(1,2)) ||  all( ch1_time_footprint > Molecule(k).Results(end,2)) )
                    
                    % Now let's define some shorthand for the necessary
                    % data from this molecule
                    
                    % Ch1 data
                    time1 = mean_decimate_array(Molecule(j).Results(bidx:end,2), decimation_factor);
                    x1 = mean_decimate_array(Molecule(j).Results(:,3), decimation_factor); y1 = mean_decimate_array(Molecule(j).Results(bidx:end,4), decimation_factor);
                    
                    % Ch2 data
                    time2 = mean_decimate_array(Molecule(k).Results(bidx:end,2), decimation_factor);
                    x2 = mean_decimate_array(Molecule(k).Results(:,3), decimation_factor); y2 = mean_decimate_array(Molecule(k).Results(bidx:end,4), decimation_factor);
                    
                    % Assume Ch1 is shorter
                    short = 1;
                    tshort = time1; tlong = time2;
                    xshort = x1; yshort = y1;
                    xlong = x2; ylong = y2;
                    if length(time1) > length(time2) % and then flip if Ch2 is shorter
                        short = 2;
                        tshort = time2; tlong = time1;
                        xlong = x1; ylong = y1;
                        xshort = x2; yshort = y2;
                    end

                        % closest beginning
                    
                    if ~isempty(tshort)

                        xdiff = [];
                        ydiff = [];
                        for a = 1:length(tshort)
                            [~,closestIndex] = min(abs(tlong-tshort(a)));
                    
                            % switch whether short is 1 or 2
                            if mod(short,2)
                                xdiff = [xdiff, xlong(closestIndex) - xshort(a)]; % x2-x1
                                ydiff = [ydiff, ylong(closestIndex) - yshort(a)]; % x2-x1
                    
                            else
                                xdiff = [xdiff, -xlong(closestIndex) + xshort(a)]; % x1-x2
                                ydiff = [ydiff, -ylong(closestIndex) + yshort(a)]; % x1-x2
                    
                            end
                        end
                        
                        startidx = length(dx_all)+1;
                        dx_all = [dx_all, xdiff];
                        dy_all = [dy_all, ydiff];
    
                        % Print a report so that we can look through them for
                        % outliers. Maybe remove automatically in the future.
                        fprintf(strcat(Molecule(j).Name, "  &  ", Molecule(k).Name , "\n"))
                        fprintf(strcat('Mean dx, dy: (', num2str(round(mean(xdiff,'omitnan'),3)) ,', ', num2str(round(mean(ydiff,'omitnan'),3)), ')\n'))
                        fprintf("\n")
                        
                        st = strfind(Molecule(j).Name, " ");
                        molnum = [molnum; str2double(Molecule(j).Name(st(end):end))];
                        molinfo(size(molinfo,1)+1,:) = [mean(xdiff,'omitnan'), mean(ydiff,'omitnan'), startidx, length(dx_all), tshort(1)];
                   end
                end

            end


        end
    end

end

% Before post-process outlier removed

dx = mean(dx_all);
dy = mean(dy_all);

stdx = std(dx_all);
stdy = std(dy_all);

% figure()
% subplot(1,2,1)
% histogram(dx_all)
% title('\Delta X before \pm 1\sigma outliers removed')
% legend(strcat("\mu: ", num2str(round(dx,2)), ", \sigma: ", num2str(round(stdx,2))))
% subplot(1,2,2)
% histogram(dy_all)
% title('\Delta Y before \pm 1\sigma outliers removed')
% legend(strcat("\mu: ", num2str(round(dy,2)), ", \sigma: ", num2str(round(stdy,2))))

figure()
subplot(1,2,1)
scatter(molinfo(:,5),molinfo(:,1), 20, 'filled')
title("X Drift")
% legend(strcat("\mu: ", num2str(round(dx,2)), ", \sigma: ", num2str(round(stdx,2))))
subplot(1,2,2)
scatter(molinfo(:,5),molinfo(:,2), 20, 'filled')
title("Y Drift")
% legend(strcat("\mu: ", num2str(round(dy,2)), ", \sigma: ", num2str(round(stdy,2))))


% After post-process outlier removed

manual_factor = round(std_max/sqrt(stdx.^2 + stdy.^2),2);
std_max_factor = min(manual_factor, std_max_factor);

removed_mols = [];
for i = size(molinfo,1):-1:1
    % use median, robust mean, instead of mean just because the noise is
    % more likely symmetric and the robust mean will be better
    if sqrt((median(dx_all) - molinfo(i,1)).^2 + (median(dy_all) - molinfo(i,2)).^2) > std_max_factor*sqrt(stdx.^2 + stdy.^2)
        dx_all(molinfo(i,3):molinfo(i,4)) = [];
        dy_all(molinfo(i,3):molinfo(i,4)) = [];
        removed_mols(size(removed_mols,1)+1,:) = molinfo(i,:);
    end
end

removed_mols

dx = mean(dx_all);
dy = mean(dy_all);

stdx = std(dx_all);
stdy = std(dy_all);


figure()
subplot(1,2,1)
scatter(molinfo(:,5),molinfo(:,1), 20, 'filled')
hold on
for i = 1:size(molnum)
    text(molinfo(i,5),molinfo(i,1)+0.25, num2str(molnum(i)))
end
if ~isempty(removed_mols)
scatter(removed_mols(:,5), removed_mols(:,1), 20, 'filled')
end
plot([0,1.1*max(molinfo(:,5))],[dx,dx],'r--')
title("X Drift w/ \pm" + num2str(std_max_factor) + " \sigma outliers removed")
ax = gca;
set(ax, 'FontName', 'Arial', 'FontSize', 10, 'TickDir', 'out', 'LineWidth', 1, 'Box', 'off', 'XColor', 'k', 'YColor', 'k');

subplot(1,2,2)
scatter(molinfo(:,5),molinfo(:,2), 20, 'filled')
hold on
for i = 1:size(molnum)
    text(molinfo(i,5),molinfo(i,2)+0.25, num2str(molnum(i)))
end
if ~isempty(removed_mols)
scatter(removed_mols(:,5), removed_mols(:,2), 20, 'filled')
end
plot([0,1.1*max(molinfo(:,5))],[dy,dy],'r--')
title("Y Drift w/ \pm" + num2str(std_max_factor) + " \sigma outliers removed")
ax = gca;
set(ax, 'FontName', 'Arial', 'FontSize', 10, 'TickDir', 'out', 'LineWidth', 1, 'Box', 'off', 'XColor', 'k', 'YColor', 'k');


% Separate option, now combined into one total figure
% figure()
% scatter(molinfo(:,1),molinfo(:,2), 20, 'filled')
% hold on
% for i = 1:size(molnum)
%     text(molinfo(i,1),molinfo(i,2)+0.25, num2str(molnum(i)))
% end
% if ~isempty(removed_mols)
% scatter(removed_mols(:,1), removed_mols(:,2), 20, 'filled')
% end
% % plot([0,1.1*max(molinfo(:,5))],[dx,dx],'r--')
% title("X-Y mean correction w/ \pm" + num2str(std_max_factor) + " \sigma outliers removed")
% xlabel('X offset correction (nm)')
% ylabel('Y offset correction (nm)')
% ax = gca;
% set(ax, 'FontName', 'Arial', 'FontSize', 10, 'TickDir', 'out', 'LineWidth', 1, 'Box', 'off', 'XColor', 'k', 'YColor', 'k');
% 
% 
% figure()
% subplot(1,2,1)
% histogram(dx_all)
% title("\Delta X after \pm" + num2str(std_max_factor) + " \sigma outliers removed")
% legend(strcat("\mu: ", num2str(round(dx,2)), ", \sigma: ", num2str(round(stdx,2))))
% ax = gca;
% set(ax, 'FontName', 'Arial', 'FontSize', 10, 'TickDir', 'out', 'LineWidth', 1, 'Box', 'off', 'XColor', 'k', 'YColor', 'k');
% subplot(1,2,2)
% histogram(dy_all)
% title("\Delta Y after \pm" + num2str(std_max_factor) + " \sigma outliers removed")
% legend(strcat("\mu: ", num2str(round(dy,2)), ", \sigma: ", num2str(round(stdy,2))))
% ax = gca;
% set(ax, 'FontName', 'Arial', 'FontSize', 10, 'TickDir', 'out', 'LineWidth', 1, 'Box', 'off', 'XColor', 'k', 'YColor', 'k');


% Figure setup
figure('Position', [200 200 600 600]);

% Define axis positions (tight layout, right-side y histogram)
scatter_pos = [0.13 0.13 0.65 0.65];  % [left bottom width height]
xhist_pos   = [0.13 0.79 0.65 0.18];  % top histogram
yhist_pos   = [0.79 0.13 0.18 0.65];  % right histogram

% Scatter plot
axScatter = axes('Position', scatter_pos);
scatter(axScatter, molinfo(:,1), molinfo(:,2), 25, 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('X offset correction (nm)');
ylabel('Y offset correction (nm)');
hold on;
% for i = 1:size(molnum)
%     text(molinfo(i,1),molinfo(i,2)+0.25, num2str(molnum(i)))
% end
if ~isempty(removed_mols)
scatter(axScatter, removed_mols(:,1), removed_mols(:,2), 20, 'filled')
end
set(axScatter, 'FontName', 'Arial', 'FontSize', 10, 'TickDir', 'out', 'LineWidth', 1, 'Box', 'off', 'XColor', 'k', 'YColor', 'k');
axis equal

% X histogram (top)
axHistX = axes('Position', xhist_pos);
histogram(axHistX, dx_all, 30, 'FaceColor', [0.5 0.5 0.5], 'BinWidth', 2);
axHistX.XAxis.Visible = 'off';
axHistX.YAxis.Visible = 'off';
axHistX.Box = 'off';
legend(strcat("\mu: ", num2str(round(dx,2)), ", \sigma: ", num2str(round(stdx,2))))
xlim(axHistX, axScatter.XLim);
set(axHistX, 'FontName', 'Arial', 'FontSize', 10, 'TickDir', 'out', 'LineWidth', 1, 'Box', 'off', 'XColor', 'k', 'YColor', 'k');

% Y histogram (right)
axHistY = axes('Position', yhist_pos);
histogram(axHistY, dy_all, 30, 'FaceColor', [0.5 0.5 0.5], 'BinWidth', 2, 'Orientation', 'horizontal');
axHistY.XAxis.Visible = 'off';
axHistY.YAxis.Visible = 'off';
axHistY.Box = 'off';
legend(strcat("\mu: ", num2str(round(dy,2)), ", \sigma: ", num2str(round(stdy,2))))
ylim(axHistY, axScatter.YLim);
set(axHistY, 'FontName', 'Arial', 'FontSize', 10, 'TickDir', 'out', 'LineWidth', 1, 'Box', 'off', 'XColor', 'k', 'YColor', 'k');

% Link axes for alignment
linkaxes([axScatter axHistX], 'x');
linkaxes([axScatter axHistY], 'y');

% Bring scatter to front
axes(axScatter);


%% Extra functions

function B = mean_decimate_array(A, decimation_factor)

% Compute the number of complete groups of decimation_factor
numGroups = floor(length(A) / decimation_factor);

% Truncate the array to only include full groups of decimation_factor
truncatedData = A(1:(numGroups * decimation_factor));

% Reshape the truncated array into 5 rows (each column corresponds to a group of decimation_factor points)
reshapedData = reshape(truncatedData, decimation_factor, []);

% Compute the mean of each column
B = mean(reshapedData, 1);