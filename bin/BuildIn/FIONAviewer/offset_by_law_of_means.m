function [don, doff] = offset_by_law_of_means(dir, fnames)
% Use law of means to just guess the offset for a set of traces. Should
% converge by the CLT, so make sure you have a sufficient number of traces.

% Parameters
% dir : (string) parent folder location of files
% fnames : (class of string(s)) filenames

% Returns
% don : (float) average difference of ch2 - ch1 in on-axis, thus correction
% to ch1
% doff : (float) average difference of ch2 - ch1 in off-axis, thus
% correction to ch1

%% We begin with the normal script for opening the data for all files and their neighbors.
% Taken verbadum from all_in_one_plot_fig_pdf

% % debugging
% fnames = {' 191838_nbh_191837_fiona.mat', ' 176368_nbh_176365_fiona.mat', ' 93072_nbh_93071_fiona.mat', ' 81725_nbh_81723_fiona.mat', ' 76838_nbh_76837_fiona.mat', ' 51677_nbh_51676_fiona.mat', ' 43311_nbh_43310_fiona.mat', ' 22703_nbh_22700_fiona.mat', ' 19056_nbh_19055_fiona.mat'};
% dir = '/Volumes/TOSHIBA EXT/MINFLUX JS/241216/241216_03/241216-110351_minflux_655_LP7p5_bg40k_555_LP4_bg15k_fies/A Neighbors/';
% 
% fnames = {' 17425_nbh_17424_fiona.mat', ' 53743_nbh_53742_fiona.mat', ' 75090_nbh_75089_fiona.mat', ' 122995_nbh_122992_fiona.mat'};
% dir = '/Volumes/TOSHIBA EXT/MINFLUX JS/241217/241217_07/241217-133308_minflux_655_LP6p5_bg40k_555_LP4_bg15k_fies/ANeighbors/';

if nargin < 1
    [fnames, dir] = uigetfile({'*.mat'},'MultiSelect','on');
end

don_all = [];
doff_all = [];

don_step_all = [];
doff_step_all = [];

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

    % find neighbor file

    jj = strfind(fname,'_fiona')-1; %nbh or initial
    j = strfind(fname,'_nbh')-1; %if there is a neighbor
    ii = 1;
    if ~isempty(j)
        ii = j+6;
    end
    
    % fname(ii:jj)
    % fname(1:j)

    if isfield(ch2data, 'time')
        time2 = ch2data.time;
        time2 = ch2data.time(~isnan(ch2data.time));
    else
        time2 = 1:size(ch2data.xy,1);
    end

    idx2 = find(~isnan(ch2data.trace(:,1)));


    % Same thing for Molecule 1

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

    end

    % Now that this file is loaded, we proceed to grab the necessary
    % information
    
    % also available in trace, but just remove the NaNs so they correspond
    % with time
    x1on = find(~isnan(ch1data.xy(:,1)));
    x2on = find(~isnan(ch2data.xy(:,1)));

    % First order mean is literally just find the average difference
    % between every nearby time point.
    % Unfortunately, there are no exact matches so we'll just have to take
    % the closest
    
    % this is grossly inefficient, but we'll do it anyway
    % assume idx1 is more populated with Nans
    short = 1;
    tshort = time1; tlong = time2;
    xshort = ch1data.xy(x1on,1); yshort = ch1data.xy(x1on,2);
    xtraceshort = ch1data.trace(x1on,3); ytraceshort = ch1data.trace_yx(x1on,3);
    xlong = ch2data.xy(x2on,1); ylong = ch2data.xy(x2on,2);
    xtracelong = ch2data.trace(x2on,3); ytracelong = ch2data.trace_yx(x2on,3);
    if length(time1) > length(time2)
        short = 2;
        tshort = time2; tlong = time1;
        xlong = ch1data.xy(x1on,1); ylong = ch1data.xy(x1on,2);
        xtracelong = ch1data.trace(x1on,3); ytracelong = ch1data.trace_yx(x1on,3);
        xshort = ch2data.xy(x2on,1); yshort = ch2data.xy(x2on,2);
        xtraceshort = ch2data.trace(x2on,3); ytraceshort = ch2data.trace_yx(x2on,3);
    end
    
    % closest beginning

    ondiff = [];
    offdiff = [];
    onstepdiff = [];
    offstepdiff = [];
    for k = 1:length(tshort)
        [~,closestIndex] = min(abs(tlong-tshort(k)));

        % switch whether short is 1 or 2
        if mod(short,2)
            ondiff = [ondiff, xlong(closestIndex) - xshort(k)]; % x2-x1
            offdiff = [offdiff, ylong(closestIndex) - yshort(k)]; % x2-x1
            
            onstepdiff = [onstepdiff, xtracelong(closestIndex) - xtraceshort(k)];
            offstepdiff = [offstepdiff, ytracelong(closestIndex) - ytraceshort(k)];

        else
            ondiff = [ondiff, -xlong(closestIndex) + xshort(k)]; % x1-x2
            offdiff = [offdiff, -ylong(closestIndex) + yshort(k)]; % x1-x2

            onstepdiff = [onstepdiff, -xtracelong(closestIndex) + xtraceshort(k)];
            offstepdiff = [offstepdiff, -ytracelong(closestIndex) + ytraceshort(k)];

        end
    end
    % figure()
    % histogram(ondiff)
    % title(fname1,"Interpreter","None")
    % figure()
    % histogram(offdiff)

    don_all = [don_all, ondiff];
    doff_all = [doff_all, offdiff];

    don_step_all = [don_step_all, onstepdiff];
    doff_step_all = [doff_step_all, offstepdiff];

    % get steps as well for future processing? Idk about this one just yet

end

sum(isnan(don_step_all))

don = mean(don_all);
doff = mean(doff_all);

stdon = std(don_all);
stdoff = std(doff_all);

figure()
histogram(don_all)
title('On-axis Difference')
legend(strcat("\mu: ", num2str(round(don,2)), ", \sigma: ", num2str(round(stdon,2))))
figure()
histogram(doff_all)
title('Off-axis Difference')
legend(strcat("\mu: ", num2str(round(doff,2)), ", \sigma: ", num2str(round(stdoff,2))))



donstep = mean(don_step_all,'omitnan');
doffstep = mean(doff_step_all);

stdon = std(don_all);
stdoff = std(doff_all);

figure()
histogram(don_all)
title('On-axis Difference')
legend(strcat("\mu: ", num2str(round(donstep,2)), ", \sigma: ", num2str(round(stdon,2))))
figure()
histogram(doff_all)
title('Off-axis Difference')
legend(strcat("\mu: ", num2str(round(doff,2)), ", \sigma: ", num2str(round(stdoff,2))))
