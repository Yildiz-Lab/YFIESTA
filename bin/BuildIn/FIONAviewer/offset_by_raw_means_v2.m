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

decimation_factor = 1;
std_max = 0.5;

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

for i = 1:length(fnames)
    try
        fname = fnames{i};
    catch
        fprintf("Finished combining xy and yx traces. Aborting... \n")
        return
    end
    
    % Now dissect filename connections

    data = load(fullfile(dir,fname));
    Molecule = data.Molecule;

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
                    time1 = mean_decimate_array(Molecule(j).Results(:,2), decimation_factor);
                    x1 = mean_decimate_array(Molecule(j).Results(:,3), decimation_factor); y1 = mean_decimate_array(Molecule(j).Results(:,4), decimation_factor);
                    
                    % Ch2 data
                    time2 = mean_decimate_array(Molecule(k).Results(:,2), decimation_factor);
                    x2 = mean_decimate_array(Molecule(k).Results(:,3), decimation_factor); y2 = mean_decimate_array(Molecule(k).Results(:,4), decimation_factor);
                    
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
                    fprintf(strcat('Mean dx, dy: (', num2str(round(mean(xdiff),3)) ,', ', num2str(round(mean(ydiff),3)), ')\n'))
                    fprintf("\n")

                    molinfo(size(molinfo,1)+1,:) = [mean(xdiff), mean(ydiff), startidx, length(dx_all)];

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

% After post-process outlier removed

removed_mols = [];
for i = size(molinfo,1):-1:1
    % use median, robust mean, instead of mean just because the noise is
    % more likely symmetric and the robust mean will be better
    if sqrt((median(dx_all) - molinfo(i,1)).^2 + (median(dy_all) - molinfo(i,2)).^2) > std_max*sqrt(stdx.^2 + stdy.^2)
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
histogram(dx_all)
title("\Delta X after \pm" + num2str(std_max) + " \sigma outliers removed")
legend(strcat("\mu: ", num2str(round(dx,2)), ", \sigma: ", num2str(round(stdx,2))))
subplot(1,2,2)
histogram(dy_all)
title("\Delta Y after \pm" + num2str(std_max) + " \sigma outliers removed")
legend(strcat("\mu: ", num2str(round(dy,2)), ", \sigma: ", num2str(round(stdy,2))))


function B = mean_decimate_array(A, decimation_factor)

% Compute the number of complete groups of decimation_factor
numGroups = floor(length(A) / decimation_factor);

% Truncate the array to only include full groups of decimation_factor
truncatedData = A(1:(numGroups * decimation_factor));

% Reshape the truncated array into 5 rows (each column corresponds to a group of decimation_factor points)
reshapedData = reshape(truncatedData, decimation_factor, []);

% Compute the mean of each column
B = mean(reshapedData, 1);