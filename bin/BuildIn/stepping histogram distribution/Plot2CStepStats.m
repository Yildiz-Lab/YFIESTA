function Plot2CStepStats(xydiff, xy_deltatxy_step)
%Parameters
% xydiff (Nx2 array) - contains all data on the separation of both channels in x and y
% xy_deltatxy_step (nx7 array) - each containing information about a changepoint step with the following format:
%   [channel, normalized trace time (s), x_head_separation (nm), y_head_separation (nm), deltat (s), deltax (nm), deltay (nm)];

% Returns:
% 1. Histograms of xy-trace difference in two color separation
% 2. Scatter plots of step size and dwell time depending on interhead separation

%% Plot options

figure()
% On-axis Inter-head separation
subplot(1,2,1)
histogram(xydiff(:,1))
title('Long-axis separation (nm)')

% Off-axis Inter-head separation
subplot(1,2,2)
histogram(xydiff(:,2))
title('Short-axis separation (nm)')


%% Next figure with line fits
x = xy_deltatxy_step(:,3); %x-head separation
y = xy_deltatxy_step(:,6); %x-step
mask = ~isnan(x .* y);
x = x(mask); y = y(mask);

% Fit a linear model (y = mx + b)
coeffs = polyfit(x, y, 1); % First-degree polynomial fit
m = coeffs(1);  % Slope
b = coeffs(2);  % Intercept

% Generate fitted y values
y_fit = polyval(coeffs, x);  % Evaluate polynomial at original x values

% Compute R-squared value
SS_total = sum((y - mean(y)).^2);  % Total sum of squares
SS_residual = sum((y - y_fit).^2);  % Residual sum of squares
R2 = 1 - (SS_residual / SS_total);  % R-squared formula

x_on = [x, y]; %assign variable for later use

figure('Position', [161 169 847 595])
% Step size on-axis
subplot(3,1,1)
scatter(xy_deltatxy_step(:,3), xy_deltatxy_step(:,6), 50, 'k', 'filled')
hold on;
plot(x, y_fit, 'r-', 'LineWidth', 2, ...
    'DisplayName', sprintf('Fit: y = %.2fx + %.2f\nR^2 = %.3f', m, b, R2)); % Line of best fit
legend('Location', 'best');  % Display legend
hold off;
ylabel('\Delta long-axis (nm)')
title('Dependence on MTBD separation')

% Step size off-axis
x = xy_deltatxy_step(:,4);
y = xy_deltatxy_step(:,7);
mask = ~isnan(x .* y);
x = x(mask); y = y(mask);

% Fit a linear model (y = mx + b)
coeffs = polyfit(x, y, 1); % First-degree polynomial fit
m = coeffs(1);  % Slope
b = coeffs(2);  % Intercept

% Generate fitted y values
y_fit = polyval(coeffs, x);  % Evaluate polynomial at original x values

% Compute R-squared value
SS_total = sum((y - mean(y)).^2);  % Total sum of squares
SS_residual = sum((y - y_fit).^2);  % Residual sum of squares
R2 = 1 - (SS_residual / SS_total);  % R-squared formula

x_off = [x, y]; %assign variable for later use

% make nature style
xlim([-64,64]); xticks(-64:8:64);
ylim([-64,64]); yticks(-64:8:64);
ax = gca;
set(ax, 'FontName', 'Arial', 'FontSize', 10, 'TickDir', 'out', 'LineWidth', 1, 'Box', 'off', 'XColor', 'k', 'YColor', 'k');

% store y-axis limits for histogram breakdown later
long_axis_range = ax.YLim;

subplot(3,1,2)
scatter(xy_deltatxy_step(:,4), xy_deltatxy_step(:,7), 50, 'k', 'filled')
hold on;
plot(x, y_fit, 'r-', 'LineWidth', 2, ...
    'DisplayName', sprintf('Fit: y = %.2fx + %.2f\nR^2 = %.3f', m, b, R2)); % Line of best fit
legend('Location', 'best');  % Display legend
hold off;
ylabel('\Delta short-axis (nm)')
title('Dependence on MTBD separation')

% make nature style
xlim([-48,48]); xticks(-48:8:48);
ylim([-48,48]); yticks(-48:8:48);
ax = gca;
set(ax, 'FontName', 'Arial', 'FontSize', 10, 'TickDir', 'out', 'LineWidth', 1, 'Box', 'off', 'XColor', 'k', 'YColor', 'k');

% store y-axis limits for histogram breakdown later
short_axis_range = ax.YLim;

% Dwell
subplot(3,1,3)
scatter(xy_deltatxy_step(:,3), xy_deltatxy_step(:,5), 50, 'k', 'filled')
ylabel('Dwell time')
xlabel('MTBD separation (nm)')

% make nature style
ax = gca;
set(ax, 'FontName', 'Arial', 'FontSize', 10, 'TickDir', 'out', 'LineWidth', 1, 'Box', 'off', 'XColor', 'k', 'YColor', 'k');



%% Other Mark DeWitt-esque plots (2012)

%% Histograms split by whether head is leading or trailing / right or left
%% Fraction of steps taken by leading vs trailing head depending on distance (binned)

% 1 is that the opposite head is reference, -1 is that current head is reference, equivalent to reflection across y-axis
% Make sure this aligns with Line 155 definition of the same variable in
% two_color_statistics.m
flip_reference_head = 1;

% 1 means that the opposite head is reference. So positive means I am the
% trailing head and negative means I am in the lead.
% 1 means that the left head is positive and the right head negative
% (opposite plot to Mark)

% This changes to Mark DeWitt if switch to -1. Make sure to do it in
% two_color_statistics.m

if flip_reference_head > 0
    trailing_idx = find(x_on(:,1) > 0);
    leading_idx = find(x_on(:,1) < 0);
    left_idx = find(x_off(:,1) > 0);
    right_idx = find(x_off(:,1) < 0);
else
    trailing_idx = find(x_on(:,1) < 0);
    leading_idx = find(x_on(:,1) > 0);
    left_idx = find(x_off(:,1) < 0);
    right_idx = find(x_off(:,1) > 0);
end


% Now let's plot the histograms

figure()
subplot(2,2,ceil(1-flip_reference_head/2))
h1 = histogram(x_on(leading_idx,2),'BinWidth',2);
title('Leading head')
ax1=gca; set(ax1,'XLim',long_axis_range);
set(ax1, 'FontName', 'Arial', 'FontSize', 10, 'TickDir', 'out', 'LineWidth', 1, 'Box', 'off', 'XColor', 'k', 'YColor', 'k');

subplot(2,2,ceil(1+flip_reference_head/2))
h2 = histogram(x_on(trailing_idx,2),'BinWidth',2);
title('Trailing head')
ax2=gca; set(ax2,'XLim',long_axis_range);
set(ax2, 'FontName', 'Arial', 'FontSize', 10, 'TickDir', 'out', 'LineWidth', 1, 'Box', 'off', 'XColor', 'k', 'YColor', 'k');

% set to the same amplitude for axis counts
set(ax1, 'YLim', [0,max([h1.Values, h2.Values])+10])
set(ax2, 'YLim', [0,max([h1.Values, h2.Values])+10])

subplot(2,2,ceil(3-flip_reference_head/2))
h3 = histogram(x_off(right_idx,2),'BinWidth',2);
title('Right head')
ax3=gca; set(ax3,'XLim',short_axis_range);
set(ax3, 'FontName', 'Arial', 'FontSize', 10, 'TickDir', 'out', 'LineWidth', 1, 'Box', 'off', 'XColor', 'k', 'YColor', 'k');

subplot(2,2,ceil(3+flip_reference_head/2))
h4 = histogram(x_off(left_idx,2),'BinWidth',2);
title('Left head')
ax4=gca; set(ax4,'XLim',short_axis_range);
set(ax4, 'FontName', 'Arial', 'FontSize', 10, 'TickDir', 'out', 'LineWidth', 1, 'Box', 'off', 'XColor', 'k', 'YColor', 'k');

% set to the same amplitude for axis counts
set(ax3, 'YLim', [0,max([h3.Values, h4.Values])+10])
set(ax4, 'YLim', [0,max([h3.Values, h4.Values])+10])


%% And now the fraction of steps taken by each population, at least in on-axis

figure('Position',[311 266 722 420])
subplot(1,2,1)
% binedges = -95:10:95; %make symmetric across 0 for ease of flipping
binedges = -60:8:60; %make symmetric across 0 for ease of flipping
last_bin_show = 5;

hold on
h_trail = histogram(x_on(trailing_idx,1),'BinEdges',binedges,'DisplayName','Trailing');
h_lead = histogram(x_on(leading_idx,1),'BinEdges',binedges,'DisplayName','Leading');
% h_trail = histogram(x_on(trailing_idx,2),'BinEdges',binedges,'DisplayName','Trailing');
% h_lead = histogram(x_on(leading_idx,2),'BinEdges',binedges,'DisplayName','Leading');
hold off
title('Population steps taken')
legend()
ax = gca;
set(ax, 'FontName', 'Arial', 'FontSize', 10, 'TickDir', 'out', 'LineWidth', 1, 'Box', 'off', 'XColor', 'k', 'YColor', 'k');

N_trail = h_trail.Values;
N_lead = h_lead.Values;

binidx = find(binedges > 0);
% binedges(binidx(1:end-1)) + (binedges(binidx(2:end)) - binedges(binidx(1:end-1)))/2
sep = binedges(binidx(1:end-1))/2 + binedges(binidx(2:end))/2;

subplot(1,2,2)
if flip_reference_head
    N_lead = fliplr(N_lead);
else
    N_trail = fliplr(N_trail);
end

% N_lead = N_lead(binidx(1:end-1));
% N_trail = N_trail(binidx(1:end-1));
N_lead = N_lead(binidx(1:last_bin_show));
N_trail = N_trail(binidx(1:last_bin_show));
sep = sep(1:last_bin_show);

bcf_trail = nan(length(N_trail),3);
bcf_lead = nan(length(N_trail),3);
for i = 1:length(N_trail)
    % beta_confidence(N_trail(i), N_lead(i))
    % if N_trail(i) + N_lead(i) > 0 % to prevent emptiness problems
        [a,b,c] = beta_confidence(N_trail(i), N_lead(i));
        bcf_trail(i,1) = a; bcf_trail(i,2) = b; bcf_trail(i,3) = c;
        [a,b,c] = beta_confidence(N_lead(i), N_trail(i));
        bcf_lead(i,1) = a; bcf_lead(i,2) = b; bcf_lead(i,3) = c;
    % end
end

hold on
errorbar(sep, bcf_trail(:,1), bcf_trail(:,2)-bcf_trail(:,1), bcf_trail(:,3)-bcf_trail(:,1), 'Color', [0 0.4470 0.7410], 'LineWidth',2, 'DisplayName','Trailing head')
errorbar(sep, bcf_lead(:,1), bcf_lead(:,2)-bcf_lead(:,1), bcf_lead(:,3)-bcf_lead(:,1), 'Color', [0.8500 0.3250 0.0980], 'LineWidth',2, 'DisplayName','Leading head')
legend('Location', 'best','AutoUpdate','off');  % Display legend

scatter(sep, bcf_trail(:,1), 50, [0 0.4470 0.7410], 'filled', 'DisplayName','')
scatter(sep, bcf_lead(:,1), 50, [0.8500 0.3250 0.0980], 'filled', 'DisplayName','')

xlabel('MTBD separation (nm)')
ylabel('Percentage of steps')
ylim([0,1])
xlim([binedges(round(end/2)+1), binedges(round(end/2)+1+last_bin_show)]);
%if including 0nm, include this:
% xlim([binedges(round(end/2)), binedges(round(end/2)+last_bin_show)]); 
xticks(sep)
ax = gca;
set(ax, 'FontName', 'Arial', 'FontSize', 10, 'TickDir', 'out', 'LineWidth', 1, 'Box', 'off', 'XColor', 'k', 'YColor', 'k');

hold off

%% Dwell times binning and averaging, adapted from file to show the dwell times vs interhead separation


% load from object associated with dependence on interhead separation

% dwellspacing - XData
% dwelltime - YData
dwellspacing = xy_deltatxy_step(:,3)';
dwelltime = xy_deltatxy_step(:,5)';

% bins = -44:8:60;
bins = -52:8:52;

dwell_times_binned = cell(1,length(bins)-1);

for j = 1:length(bins)-1 % bin your dwells
    mask1 = and(dwellspacing > bins(j), dwellspacing < bins(j+1));
    dwell_times_binned{j} = dwelltime(mask1)*1000; %convert to ms
end



try

figure()
hold on

mean_dwell_time = nan(length(bins)-1,3);



for j = 1:length(bins)-1 % now fit distribution and get the tau value with confidence intervals
    
    pd = fitdist(dwell_times_binned{j}', 'Exponential');
    % optional, remove outliers 2025/08/13

    % plot original for comparison
    % scatter((bins(j+1)+bins(j))/2 * ones(length(dwell_times_binned{j}),1), dwell_times_binned{j}, 20, 'r', 'filled')

    % Quantile cutoff
    % cutoff = quantile(dwell_times_binned{j}, 0.9);
    % dwell_times_binned{j} = dwell_times_binned{j}(dwell_times_binned{j} <= cutoff);
    % pd = fitdist(dwell_times_binned{j}', 'Exponential');

    % Remove based off pd estimate
    % dwell_times_binned{j}(dwell_times_binned{j} > 5*pd.mu) = [];
    
    % Remove using a double exponential fit
    % Also uncomment function at the end of file
    % Example data (replace with dwell_times_binned{j})
    % t = exprnd(1, 200, 1);               % fast process
    % t = [t; exprnd(5, 50, 1)];           % plus slow process
    % t = dwell_times_binned{j};
    % 
    % % Initial guesses: [A, tau1, tau2]
    % initParams = [0.1, 20, 100];
    % 
    % % Likelihood function for mixture of two exponentials
    % negLogLik = @(params) -sum( log( ...
    %     params(1) .* (1/params(2)) .* exp(-t/params(2)) + ...
    %     (1 - params(1)) .* (1/params(3)) .* exp(-t/params(3)) ...
    % ) );
    % 
    % % Constrain A between 0 and 1, taus > 0
    % opts = optimset('fminsearch');
    % % opts.Display = 'iter';
    % params_hat = fminsearch(@(p) penalizedNegLL(p, negLogLik), initParams, opts);
    % 
    % A_est     = params_hat(1);
    % tau1_est  = params_hat(2);
    % tau2_est  = params_hat(3);
    % 
    % fprintf('A = %.3f, tau1 = %.3f, tau2 = %.3f\n', A_est, tau1_est, tau2_est);
    
    scatter((bins(j+1)+bins(j))/2 * ones(length(dwell_times_binned{j}),1), dwell_times_binned{j}, 20, 'k', 'filled')
    
    ci = paramci(pd);
    mean_dwell_time(j,1) = pd.mu; mean_dwell_time(j,2:3) = ci';
end

scatter((bins(1:end-1)+bins(2:end))/2, mean_dwell_time(:,1), 20, 'r', 'filled')
errorbar((bins(1:end-1)+bins(2:end))/2, mean_dwell_time(:,1), mean_dwell_time(:,1) - mean_dwell_time(:,2), -mean_dwell_time(:,1) + mean_dwell_time(:,3), 'r', 'LineStyle', 'none')


ax = gca;
ax.XLim = [min(bins), max(bins)];
xlabel("Interhead separation (nm)");
ylabel("Dwell Time (ms)");
set(ax, ...
        'FontName', 'Arial', ...
        'FontSize', 10, ...
        'TickDir', 'out', ...
        'LineWidth', 1, ...
        'Box', 'off', ...
        'XColor', 'k', ...
        'YColor', 'k');
catch
    fprintf('Skipping rate plotting due to insufficient data (Check Line 278-354 in Plot2CStepStats.m \n')
end

figure()
hold on
[val0, idx0] = min(abs((bins(1:end-1)+bins(2:end))/2));
if val0 < 0
    idx0 = idx0-1;
end
% Leading
scatter(-(bins(1:idx0)+bins(2:idx0+1))/2, mean_dwell_time(1:idx0,1), 40, [0.8500 0.3250 0.0980], 'filled', 'DisplayName', '')
errorbar(-(bins(1:idx0)+bins(2:idx0+1))/2, mean_dwell_time(1:idx0,1), mean_dwell_time(1:idx0,1) - mean_dwell_time(1:idx0,2), -mean_dwell_time(1:idx0,1) + mean_dwell_time(1:idx0,3), 'Color', [0.8500 0.3250 0.0980], 'LineStyle', 'none', 'DisplayName', 'Leading head', 'LineWidth',2)
% Trailing
scatter((bins(idx0+1:end-1)+bins(idx0+2:end))/2, mean_dwell_time(idx0+1:end,1), 40, [0 0.4470 0.7410], 'filled', 'DisplayName', '')
errorbar((bins(idx0+1:end-1)+bins(idx0+2:end))/2, mean_dwell_time(idx0+1:end,1), mean_dwell_time(idx0+1:end,1) - mean_dwell_time(idx0+1:end,2), -mean_dwell_time(idx0+1:end,1) + mean_dwell_time(idx0+1:end,3), 'Color', [0 0.4470 0.7410], 'LineStyle', 'none', 'DisplayName', 'Trailing head', 'LineWidth',2)

xticks((bins(idx0:end-1)+bins(idx0+1:end))/2)
yticks(5:5:25)

legend()

ax = gca;
ax.XLim = [val0 - (bins(2)-bins(1))/2, max(bins)];
xlabel("Interhead separation (nm)");
ylabel("Dwell Time (ms)");
set(ax, ...
        'FontName', 'Arial', ...
        'FontSize', 10, ...
        'TickDir', 'out', ...
        'LineWidth', 1, ...
        'Box', 'off', ...
        'XColor', 'k', ...
        'YColor', 'k');


%% And now the fraction of steps taken by each population in off-axis

figure('Position',[311 266 722 420])
subplot(1,2,1)
% binedges = -95:10:95; %make symmetric across 0 for ease of flipping
binedges = -32.5:5:32.5; %make symmetric across 0 for ease of flipping
last_bin_show = 7;

hold on
% h_trail = histogram(x_off(left_idx,2),'BinEdges',binedges,'DisplayName','Left');
% h_lead = histogram(x_off(right_idx,2),'BinEdges',binedges,'DisplayName','Right');
h_trail = histogram(x_off(left_idx,1),'BinEdges',binedges,'DisplayName','Left');
h_lead = histogram(x_off(right_idx,1),'BinEdges',binedges,'DisplayName','Right');

hold off
title('Population steps taken')
legend()
ax = gca;
set(ax, 'FontName', 'Arial', 'FontSize', 10, 'TickDir', 'out', 'LineWidth', 1, 'Box', 'off', 'XColor', 'k', 'YColor', 'k');

N_trail = h_trail.Values;
N_lead = h_lead.Values;

binidx = find(binedges > -3);
% binedges(binidx(1:end-1)) + (binedges(binidx(2:end)) - binedges(binidx(1:end-1)))/2
sep = binedges(binidx(1:end-1))/2 + binedges(binidx(2:end))/2;

subplot(1,2,2)
if flip_reference_head
    N_lead = fliplr(N_lead);
else
    N_trail = fliplr(N_trail);
end

% N_lead = N_lead(binidx(1:end-1));
% N_trail = N_trail(binidx(1:end-1));
N_lead = N_lead(binidx(1:last_bin_show));
N_trail = N_trail(binidx(1:last_bin_show));
sep = sep(1:last_bin_show);

bcf_trail = nan(length(N_trail),3);
bcf_lead = nan(length(N_trail),3);
for i = 1:length(N_trail)
    % beta_confidence(N_trail(i), N_lead(i))
    % if N_trail(i) + N_lead(i) > 0 % to prevent emptiness problems
        [a,b,c] = beta_confidence(N_trail(i), N_lead(i));
        bcf_trail(i,1) = a; bcf_trail(i,2) = b; bcf_trail(i,3) = c;
        [a,b,c] = beta_confidence(N_lead(i), N_trail(i));
        bcf_lead(i,1) = a; bcf_lead(i,2) = b; bcf_lead(i,3) = c;
    % end
end

hold on
errorbar(sep, bcf_trail(:,1), bcf_trail(:,2)-bcf_trail(:,1), bcf_trail(:,3)-bcf_trail(:,1), 'Color', [0 0.4470 0.7410], 'LineWidth',2, 'DisplayName','Left MTBD')
errorbar(sep, bcf_lead(:,1), bcf_lead(:,2)-bcf_lead(:,1), bcf_lead(:,3)-bcf_lead(:,1), 'Color', [0.8500 0.3250 0.0980], 'LineWidth',2, 'DisplayName','Right MTBD')
legend('Location', 'best','AutoUpdate','off');  % Display legend

scatter(sep, bcf_trail(:,1), 50, [0 0.4470 0.7410], 'filled', 'DisplayName','')
scatter(sep, bcf_lead(:,1), 50, [0.8500 0.3250 0.0980], 'filled', 'DisplayName','')

xlabel('MTBD separation (nm)')
ylabel('Percentage of steps')
ylim([0,1])
% xlim([binedges(round(end/2)+1), binedges(round(end/2)+1+last_bin_show)]);
%if including 0nm, include this:
xlim([binedges(round(end/2)), binedges(round(end/2)+last_bin_show)]); 
xticks(sep)
ax = gca;
set(ax, 'FontName', 'Arial', 'FontSize', 10, 'TickDir', 'out', 'LineWidth', 1, 'Box', 'off', 'XColor', 'k', 'YColor', 'k');

hold off


%% Dwell times binning and averaging, adapted from file to show the dwell times vs interhead separation


% load from object associated with dependence on interhead separation

% dwellspacing - XData
% dwelltime - YData
dwellspacing = xy_deltatxy_step(:,4)';
dwelltime = xy_deltatxy_step(:,5)';

% bins = -44:8:60;
bins = -32.5:5:32.5;

dwell_times_binned = cell(1,length(bins)-1);

for j = 1:length(bins)-1 % bin your dwells
    mask1 = and(dwellspacing > bins(j), dwellspacing < bins(j+1));
    dwell_times_binned{j} = dwelltime(mask1)*1000; %convert to ms
end

figure()

try

    subplot(2,1,1)
    hold on
mean_dwell_time = nan(length(bins)-1,3);

for j = 1:length(bins)-1 % now fit distribution and get the tau value with confidence intervals
    
    pd = fitdist(dwell_times_binned{j}', 'Exponential');
    % optional, remove outliers 2025/08/13

    % plot original for comparison
    % scatter((bins(j+1)+bins(j))/2 * ones(length(dwell_times_binned{j}),1), dwell_times_binned{j}, 20, 'r', 'filled')

    % Quantile cutoff
    % cutoff = quantile(dwell_times_binned{j}, 0.9);
    % dwell_times_binned{j} = dwell_times_binned{j}(dwell_times_binned{j} <= cutoff);
    % pd = fitdist(dwell_times_binned{j}', 'Exponential');

    % Remove based off pd estimate
    % dwell_times_binned{j}(dwell_times_binned{j} > 5*pd.mu) = [];
    
    % Remove using a double exponential fit
    % Also uncomment function at the end of file
    % Example data (replace with dwell_times_binned{j})
    % t = exprnd(1, 200, 1);               % fast process
    % t = [t; exprnd(5, 50, 1)];           % plus slow process
    % t = dwell_times_binned{j};
    % 
    % % Initial guesses: [A, tau1, tau2]
    % initParams = [0.1, 20, 100];
    % 
    % % Likelihood function for mixture of two exponentials
    % negLogLik = @(params) -sum( log( ...
    %     params(1) .* (1/params(2)) .* exp(-t/params(2)) + ...
    %     (1 - params(1)) .* (1/params(3)) .* exp(-t/params(3)) ...
    % ) );
    % 
    % % Constrain A between 0 and 1, taus > 0
    % opts = optimset('fminsearch');
    % % opts.Display = 'iter';
    % params_hat = fminsearch(@(p) penalizedNegLL(p, negLogLik), initParams, opts);
    % 
    % A_est     = params_hat(1);
    % tau1_est  = params_hat(2);
    % tau2_est  = params_hat(3);
    % 
    % fprintf('A = %.3f, tau1 = %.3f, tau2 = %.3f\n', A_est, tau1_est, tau2_est);
    
    scatter((bins(j+1)+bins(j))/2 * ones(length(dwell_times_binned{j}),1), dwell_times_binned{j}, 20, 'k', 'filled')
    
    ci = paramci(pd);
    mean_dwell_time(j,1) = pd.mu; mean_dwell_time(j,2:3) = ci';
end

scatter((bins(1:end-1)+bins(2:end))/2, mean_dwell_time(:,1), 20, 'r', 'filled')
errorbar((bins(1:end-1)+bins(2:end))/2, mean_dwell_time(:,1), mean_dwell_time(:,1) - mean_dwell_time(:,2), -mean_dwell_time(:,1) + mean_dwell_time(:,3), 'r', 'LineStyle', 'none')


ax = gca;
ax.XLim = [min(bins), max(bins)];
xlabel("Interhead separation (nm)");
ylabel("Dwell Time (ms)");
set(ax, ...
        'FontName', 'Arial', ...
        'FontSize', 10, ...
        'TickDir', 'out', ...
        'LineWidth', 1, ...
        'Box', 'off', ...
        'XColor', 'k', ...
        'YColor', 'k');
catch
    fprintf('Skipping rate plotting due to insufficient data (Check Line 278-354 in Plot2CStepStats.m \n')
end

subplot(2,1,2)
hold on
[val0, idx0] = min(abs((bins(1:end-1)+bins(2:end))/2));
if val0 < 0
    idx0 = idx0-1;
end
% Right
scatter(-(bins(1:idx0)+bins(2:idx0+1))/2, mean_dwell_time(1:idx0,1), 40, [0.8500 0.3250 0.0980], 'filled', 'DisplayName', '')
errorbar(-(bins(1:idx0)+bins(2:idx0+1))/2, mean_dwell_time(1:idx0,1), mean_dwell_time(1:idx0,1) - mean_dwell_time(1:idx0,2), -mean_dwell_time(1:idx0,1) + mean_dwell_time(1:idx0,3), 'Color', [0.8500 0.3250 0.0980], 'LineStyle', 'none', 'DisplayName', 'Right head', 'LineWidth',2)
% Left
scatter((bins(idx0+1:end-1)+bins(idx0+2:end))/2, mean_dwell_time(idx0+1:end,1), 40, [0 0.4470 0.7410], 'filled', 'DisplayName', '')
errorbar((bins(idx0+1:end-1)+bins(idx0+2:end))/2, mean_dwell_time(idx0+1:end,1), mean_dwell_time(idx0+1:end,1) - mean_dwell_time(idx0+1:end,2), -mean_dwell_time(idx0+1:end,1) + mean_dwell_time(idx0+1:end,3), 'Color', [0 0.4470 0.7410], 'LineStyle', 'none', 'DisplayName', 'Left head', 'LineWidth',2)

xticks((bins(idx0:end-1)+bins(idx0+1:end))/2)
yticks(5:5:25)

legend()

ax = gca;
ax.XLim = [val0 - (bins(2)-bins(1))/2, max(bins)];
xlabel("Interhead separation (nm)");
ylabel("Dwell Time (ms)");
set(ax, ...
        'FontName', 'Arial', ...
        'FontSize', 10, ...
        'TickDir', 'out', ...
        'LineWidth', 1, ...
        'Box', 'off', ...
        'XColor', 'k', ...
        'YColor', 'k');

%% Determine if there is one head that takes more steps (has higher overlal stepping rate)
% Is there a directional dependence on this.

% Want a correlation plot where we ask: If I am the left or right head, am I
% trailing or leading and vice versa. Maybe just a heat map. Rather than
% quadrants, will do nonants to set a threshold for distances.

boxsize = 12.5; %nm
correlation_array = zeros(4,4,2); %First row percentage, second row mean rate


% need to use the original data
mask = ~isnan(xy_deltatxy_step(:,3) .* xy_deltatxy_step(:,4));
xy_deltatxy_step_filtered = xy_deltatxy_step(mask,:); % just get rid of any steps with nan's
idx_onsteps = ~isnan(xy_deltatxy_step_filtered(:,6));
idx_offsteps = ~isnan(xy_deltatxy_step_filtered(:,7));

figure()
s = sqrt(xy_deltatxy_step_filtered(:,3).^2 + xy_deltatxy_step_filtered(:,4).^2);
histogram(s, 'DisplayName', sprintf('mu = %.2f \n sigma = %.2f', mean(s), std(s)))
legend()

ax = gca;
% ax.XLim = [val0 - (bins(2)-bins(1))/2, max(bins)];
xlabel("MTBD separation (nm)");
ylabel("Counts");
title("MTBD Separation")
set(ax, ...
        'FontName', 'Arial', ...
        'FontSize', 10, ...
        'TickDir', 'out', ...
        'LineWidth', 1, ...
        'Box', 'off', ...
        'XColor', 'k', ...
        'YColor', 'k');



if flip_reference_head > 0
    trailing_idx = find(xy_deltatxy_step_filtered(:,3) > boxsize);
    leading_idx = find(xy_deltatxy_step_filtered(:,3) < -boxsize);
    trail_within_box_idx = find(and(xy_deltatxy_step_filtered(:,3) > 0, xy_deltatxy_step_filtered(:,3) < boxsize));
    lead_within_box_idx = find(and(xy_deltatxy_step_filtered(:,3) > -boxsize, xy_deltatxy_step_filtered(:,3) < 0));
    left_idx = find(xy_deltatxy_step_filtered(:,4) > boxsize);
    right_idx = find(xy_deltatxy_step_filtered(:,4) < -boxsize);
    left_within_box_idx = find(and(xy_deltatxy_step_filtered(:,4) > 0, xy_deltatxy_step_filtered(:,4) < boxsize));
    right_within_box_idx = find(and(xy_deltatxy_step_filtered(:,4) > -boxsize, xy_deltatxy_step_filtered(:,4) < 0));
else
    trailing_idx = find(xy_deltatxy_step_filtered(:,3) < boxsize);
    leading_idx = find(xy_deltatxy_step_filtered(:,3) > -boxsize);
    lead_within_box_idx = find(and(xy_deltatxy_step_filtered(:,3) > 0, xy_deltatxy_step_filtered(:,3) < boxsize));
    trail_within_box_idx = find(and(xy_deltatxy_step_filtered(:,3) > -boxsize, xy_deltatxy_step_filtered(:,3) < 0));
    left_idx = find(xy_deltatxy_step_filtered(:,4) < boxsize);
    right_idx = find(xy_deltatxy_step_filtered(:,4) > -boxsize);
    right_within_box_idx = find(and(xy_deltatxy_step_filtered(:,4) > 0, xy_deltatxy_step_filtered(:,4) < boxsize));
    left_within_box_idx = find(and(xy_deltatxy_step_filtered(:,4) > -boxsize, xy_deltatxy_step_filtered(:,4) < 0));
end
% on_within_box_idx = find(and(xy_deltatxy_step_filtered(:,3) > -boxsize, xy_deltatxy_step_filtered(:,3) < boxsize));
% off_within_box_idx = find(and(xy_deltatxy_step_filtered(:,4) > -boxsize, xy_deltatxy_step_filtered(:,4) < boxsize));

% maybe want to separate whether the steps are on or off axis... but for
% now eh whatever.

% correlation array goes from leading left to right and so on
correlation_array(1,1,1) = length(intersect(leading_idx, left_idx));
correlation_array(1,4,1) = length(intersect(leading_idx, right_idx));
correlation_array(4,1,1) = length(intersect(trailing_idx, left_idx));
correlation_array(4,4,1) = length(intersect(trailing_idx, right_idx));

% Now within the near box
% correlation_array(1,2,1) = length(intersect(leading_idx, off_within_box_idx));
% correlation_array(2,1,1) = length(intersect(on_within_box_idx, left_idx));
% correlation_array(2,4,1) = length(intersect(on_within_box_idx, right_idx));
% correlation_array(3,2,1) = length(intersect(trailing_idx, off_within_box_idx));
% and completely unsure
% correlation_array(2,2,1) = length(intersect(on_within_box_idx, off_within_box_idx));

correlation_array(1,2,1) = length(intersect(leading_idx, left_within_box_idx));
correlation_array(1,3,1) = length(intersect(leading_idx, right_within_box_idx));
correlation_array(4,2,1) = length(intersect(trailing_idx, left_within_box_idx));
correlation_array(4,3,1) = length(intersect(trailing_idx, right_within_box_idx));

correlation_array(2,1,1) = length(intersect(lead_within_box_idx, left_idx));
correlation_array(3,1,1) = length(intersect(trail_within_box_idx, left_idx));
correlation_array(2,4,1) = length(intersect(lead_within_box_idx, right_idx));
correlation_array(3,4,1) = length(intersect(trail_within_box_idx, right_idx));

correlation_array(2,2,1) = length(intersect(lead_within_box_idx, left_within_box_idx));
correlation_array(2,3,1) = length(intersect(lead_within_box_idx, right_within_box_idx));
correlation_array(3,2,1) = length(intersect(trail_within_box_idx, left_within_box_idx));
correlation_array(3,3,1) = length(intersect(trail_within_box_idx, right_within_box_idx));

% Now for the rates (for now just do mean, we can do actual exp decay if we see
% a difference)
correlation_array(1,1,2) = mean(xy_deltatxy_step_filtered(intersect(leading_idx, left_idx),5));
correlation_array(1,4,2) = mean(xy_deltatxy_step_filtered(intersect(leading_idx, right_idx),5));
correlation_array(4,1,2) = mean(xy_deltatxy_step_filtered(intersect(trailing_idx, left_idx),5));
correlation_array(4,4,2) = mean(xy_deltatxy_step_filtered(intersect(trailing_idx, right_idx),5));

% Now within the near box
% correlation_array(1,2,2) = mean(xy_deltatxy_step_filtered(intersect(leading_idx, off_within_box_idx),5));
% correlation_array(2,1,2) = mean(xy_deltatxy_step_filtered(intersect(on_within_box_idx, left_idx),5));
% correlation_array(2,4,2) = mean(xy_deltatxy_step_filtered(intersect(on_within_box_idx, right_idx),5));
% correlation_array(3,2,2) = mean(xy_deltatxy_step_filtered(intersect(trailing_idx, off_within_box_idx),5));
% and completely unsure
% correlation_array(2,2,2) = mean(xy_deltatxy_step_filtered(intersect(on_within_box_idx, off_within_box_idx),5));

correlation_array(1,2,2) = mean(xy_deltatxy_step_filtered(intersect(leading_idx, left_within_box_idx),5));
correlation_array(1,3,2) = mean(xy_deltatxy_step_filtered(intersect(leading_idx, right_within_box_idx),5));
correlation_array(4,2,2) = mean(xy_deltatxy_step_filtered(intersect(trailing_idx, left_within_box_idx),5));
correlation_array(4,3,2) = mean(xy_deltatxy_step_filtered(intersect(trailing_idx, right_within_box_idx),5));

correlation_array(2,1,2) = mean(xy_deltatxy_step_filtered(intersect(lead_within_box_idx, left_idx),5));
correlation_array(3,1,2) = mean(xy_deltatxy_step_filtered(intersect(trail_within_box_idx, left_idx),5));
correlation_array(2,4,2) = mean(xy_deltatxy_step_filtered(intersect(lead_within_box_idx, right_idx),5));
correlation_array(3,4,2) = mean(xy_deltatxy_step_filtered(intersect(trail_within_box_idx, right_idx),5));

correlation_array(2,2,2) = mean(xy_deltatxy_step_filtered(intersect(lead_within_box_idx, left_within_box_idx),5));
correlation_array(2,3,2) = mean(xy_deltatxy_step_filtered(intersect(lead_within_box_idx, right_within_box_idx),5));
correlation_array(3,2,2) = mean(xy_deltatxy_step_filtered(intersect(trail_within_box_idx, left_within_box_idx),5));
correlation_array(3,3,2) = mean(xy_deltatxy_step_filtered(intersect(trail_within_box_idx, right_within_box_idx),5));

correlation_array
% close all
%% Autocorrelations in time
% find points to split molecules
molidx = find(xy_deltatxy_step_filtered(2:end,2) - xy_deltatxy_step_filtered(1:end-1,2) < 0); % if we have a negative time, we are at a new ref molecule
molidx = [1; molidx; size(xy_deltatxy_step_filtered,1)];

% ugh a triple for loop oof
autocorr = cell(1,30);
ch1sum = []; ch2sum = []; chtotsum = []; time_diff = [];
sliding_window = [];

r = []; theta = [];
for j = 1:length(molidx)-1
    deltatxystep = xy_deltatxy_step_filtered(molidx(j):molidx(j+1),:);
    
    % before anything, make r, theta for separation of ch1 to ch2
    % deltatxystep has column 3 and 4 containing on, off axis separation
    % appropriately. Let's just calculate it at every step.
    r = [r; sqrt(deltatxystep(:,3).^2 + deltatxystep(:,4).^2)];
    theta_add = atan(deltatxystep(:,4)./deltatxystep(:,3)) - pi*floor((sign(deltatxystep(:,3))-1)/2);
    theta = [theta; theta_add];

    % To start we should move everything into the ch1 reference
    % This is just the xsep, ysep, not the steps
    ch2mask = find(deltatxystep(:,1) == 2);
    % deltatxystep(ch2mask,[3:4,6:7]) = -deltatxystep(ch2mask,[3:4,6:7]);
    deltatxystep(ch2mask,3:4) = -deltatxystep(ch2mask,3:4);
    
    

    % could also add a buffer for separation. i.e. if the sep is too small,
    % say 5nm, then let's just say it is the same sign as the previous.
    % This isn't doing anything though...weirdly
    buffer_dist = 5;
    deltatxystep = [deltatxystep(1,:); deltatxystep]; %add first row for edge case
    too_small_idx = find(abs(deltatxystep(2:end,4)) < buffer_dist); too_small_idx = too_small_idx+1;
    for c = 1:length(too_small_idx) %have to do a for loop for iterative purposes
        deltatxystep(too_small_idx(c),4) = sign(deltatxystep(too_small_idx(c)-1,4));
    end
    too_small_idx = find(abs(deltatxystep(2:end,3)) < buffer_dist); too_small_idx = too_small_idx+1;
    for c = 1:length(too_small_idx) %have to do a for loop for iterative purposes
        deltatxystep(too_small_idx(c),3) = sign(deltatxystep(too_small_idx(c)-1,3));
    end
    too_small_idx = find(abs(deltatxystep(2:end,6)) < buffer_dist); too_small_idx = too_small_idx+1;
    for c = 1:length(too_small_idx) %have to do a for loop for iterative purposes
        deltatxystep(too_small_idx(c),6) = sign(deltatxystep(too_small_idx(c)-1,6));
    end
    too_small_idx = find(abs(deltatxystep(2:end,7)) < buffer_dist); too_small_idx = too_small_idx+1;
    for c = 1:length(too_small_idx) %have to do a for loop for iterative purposes
        deltatxystep(too_small_idx(c),7) = sign(deltatxystep(too_small_idx(c)-1,7));
    end
    deltatxystep(1,:) = []; %delete first row
    
    % now do sliding correlation by if they are close enough to each other.
    % In this case we will just do signs
    % G(1,1) = G(-1,-1) = 1;  G(1,-1) = G(-1,1) = 0
    for k = 0:length(autocorr)-1
        % first channel that is the stepper - 
        channel = -abs( deltatxystep(1:end-k,1) - deltatxystep(1+k:end,1) ) + 1;
        % then we do the other channels just based off signs autocorr
        pos_and_step = -abs( sign(deltatxystep(1:end-k,[3:4,6:7])) - sign(deltatxystep(1+k:end,[3:4,6:7])) )/2 + 1;
        temp = autocorr{k+1};
        autocorr{k+1} = [temp; channel, pos_and_step];
        
        if k == 1
            % deltatxystep
            x = [channel, pos_and_step];
            start_end_idx = find_run_regions_idx(x(:,3), 1, 5);
            % now we ask for information about time and distance during this run,
            % through the reference frame knowing the sum of all x and y steps tell
            % us the translocation of both motors
            
            ch1idx = find(deltatxystep(:,1) == 1); ch2idx = find(deltatxystep(:,1) == 2);
            
            % Now that we have this filter, we will ask for the sum of all columns in
            % between the start and the end idx separated by channels
        
            for m = 1:size(start_end_idx,1)
                mask = start_end_idx(m,1):start_end_idx(m,2)-1;
                % we subtract one because the last step is where the cross happens, and so the run has ended
                
                % calculate delta - initial condition in both channels
                mask_ch1idx = intersect(mask,ch1idx); mask_ch2idx = intersect(mask,ch2idx);
                fprintf(strcat("Number ch1: ", num2str(length(mask_ch1idx)), " Number ch2: ", num2str(length(mask_ch2idx)), "\n"))
                
                ch1sum = [ch1sum; sum(deltatxystep(mask_ch1idx,:),1), mean(x(mask_ch1idx,:),1)]; %- deltatxystep(mask_ch1idx(1),3:4);
                ch2sum = [ch2sum; sum(deltatxystep(mask_ch2idx,:),1), mean(x(mask_ch2idx,:),1)]; %- deltatxystep(mask_ch2idx(1),3:4);
                chtotsum = [chtotsum; sum(deltatxystep(mask,:),1), mean(x(mask,:),1)];
                time_diff = [time_diff; deltatxystep(mask(end),2) - deltatxystep(mask(1),2)];
                
            end
            
            % This tells us in the stretches, what is the position of the
            % relative heads
            % sign(ch1sum(end-size(start_end_idx,1)+1:end,:))
            % sign(ch2sum(end-size(start_end_idx,1)+1:end,:))
            sign(chtotsum(end-size(start_end_idx,1)+1:end,:))
            % ch1sum(end-size(start_end_idx,1)+1:end,:)
            % ch2sum(end-size(start_end_idx,1)+1:end,:)

            % Let's report some overall stats.
            % First the short axis ratio, which we have already done "run
            % stuff" based on short axis so we don't have to sum this
            % fprintf(strcat("short-axis Ratio: ", num2str( sum(sign(ch1sum(end-size(start_end_idx,1)+1:end,4)) < 1) / size(start_end_idx,1) ), '\n') )
            fprintf(strcat("short-axis Ratio: ", num2str( sum(sign(chtotsum(end-size(start_end_idx,1)+1:end,4)) < 1) / size(start_end_idx,1) ), '\n') )
            
            % Then the long-axis which depends 
            % fprintf(strcat("long-axis Ratio: ", num2str( ( sum(sign(ch1sum(end-size(start_end_idx,1)+1:end,3)) < 1) + sum(sign(ch2sum(end-size(start_end_idx,1)+1:end,3)) < 1) ) / 2 / size(start_end_idx,1) ), '\n') )
            fprintf(strcat("long-axis Ratio: ", num2str( sum(sign(chtotsum(end-size(start_end_idx,1)+1:end,3)) < 1) / size(start_end_idx,1) ), '\n') )

            % fprintf(strcat("short-axis * long-axis corr: ", num2str( ( sum( sign(ch1sum(end-size(start_end_idx,1)+1:end,3)) .* sign(ch1sum(end-size(start_end_idx,1)+1:end,4)) > 0) + sum( sign(ch2sum(end-size(start_end_idx,1)+1:end,3)) .* sign(ch2sum(end-size(start_end_idx,1)+1:end,4)) > 0) ) / 2 / size(start_end_idx,1) ), '\n') )
            fprintf(strcat("short-axis * long-axis corr: ", num2str( sum( sign(chtotsum(end-size(start_end_idx,1)+1:end,3)) .* sign(chtotsum(end-size(start_end_idx,1)+1:end,4)) > 0) / size(start_end_idx,1) ), '\n') )
            fprintf('\n')
        end

        % %autocorrelation for window of size below
        % if k == 1
        %     y = [channel, pos_and_step]; %autocorr{k+1}; %uncomment this for entire array and do it outside the loop
        %     window_size = 1;
        %     storage = nan(size(y,1)-window_size, size(deltatxystep,2)+size(y,2));
        %     for s = 1:(size(y,1)-window_size)
        %         storage(s,1:7) = mean(deltatxystep(s:s+window_size,:),'omitnan');
        %         storage(s,8:12) = mean(y(s:s+window_size-1,:),'omitnan');
        %     end
        %     % Now idea is find segments where there is "lots of crossings"
        %     % vs not and in the segments where there isn't confusion what
        %     % is dominating?
        % 
        %     not_crossing_much = find(storage(:,10) > 0.5);
        %     crossing_much = find(storage(:,10) <= 0.5);
        %     % length(not_crossing_much)
        %     % length(crossing_much)
        %     fprintf(strcat('runs w/o cross: ', num2str(mean(sign(storage(not_crossing_much,4)))), '\n'));
        %     fprintf(strcat('crosses: ', num2str(mean(sign(storage(crossing_much,4)))), '\n'));
        % 
        % end

    end

end

% r
% theta
plot_polar_conversion(r, theta, 24)

v1 = ch1sum(:,6)./time_diff;
v2 = ch2sum(:,6)./time_diff;

figure()
subplot(1,2,1)
hold on
histogram(v1,'BinWidth',25)
histogram(v2,'BinWidth',25)

subplot(1,2,2)
histogram((v1+v2)/2,'BinWidth',25)

% now let's just make a plot for average "autocorrelation" in every asset
% fixed
% deltatxystep(:,[1:4,6:7])

figure()
hold on
labels = {"channel", "xsep", "ysep", "xstep", "ystep"};
for l = 1:size(autocorr{1},2)
    g = nan(1,length(autocorr));
    for m = 1:length(g)
        x = autocorr{m};
        g(m) = mean(x(:,l),'omitnan');
    end
    
    try
        x = autocorr{2};
        x = x(:,l);
        tau = fit_exp_run_lengths(x(~isnan(x)), 1);
        fprintf(strcat(labels{l}, " tau (steps) : ", num2str(round(1/tau,2)), "\n"))
        plot(0:length(autocorr)-1,g,'DisplayName', strcat(labels{l}, " \tau : ", num2str(round(1/tau,2)), " (steps)"))
    catch
        plot(0:length(autocorr)-1,g,'DisplayName',labels{l})
    end
end
legend()
ax = gca();
set(ax, 'FontName', 'Arial', 'FontSize', 10, 'TickDir', 'out', ...
        'LineWidth', 1, 'Box', 'off', 'XColor', 'k', 'YColor', 'k');





% 
% % --- Helper function to keep parameters in bounds ---
% function nll = penalizedNegLL(p, nllFunc)
%     if p(1) < 0 || p(1) > 1 || p(2) <= 0 || p(3) <= 0
%         nll = Inf;
%     else
%         nll = nllFunc(p);
%     end

function [lambda_hat, run_lengths] = fit_exp_run_lengths(binary_array, type)
    % fit_exp_run_lengths: Fit an exponential distribution to lengths of
    % consecutive type
    % Input:
    %   binary_array - Vector of 0s and 1s
    %   type - 0 or 1 depending on if you want to know the runs of zeros or
    %   ones
    % Output:
    %   lambda_hat   - Estimated rate parameter of exponential distribution
    %   run_lengths  - Vector of run lengths of ones

    % Ensure input is a row vector and pad binary array
    binary_array = [mod(type+1,2) binary_array(:)' mod(type+1,2)];

    % Find transitions
    diff_array = diff(binary_array);
    if type
        start_idx = find(diff_array == 1);  % Start of runs
        end_idx = find(diff_array == -1); % End of runs
    else
        start_idx = find(diff_array == -1);  % Start of runs
        end_idx = find(diff_array == 1); % End of runs
    end
    

    % Run lengths
    run_lengths = end_idx - start_idx;

    if isempty(run_lengths)
        lambda_hat = NaN;
        warning('No runs of ones found.');
        return;
    end

    % Fit exponential using MATLAB's fitdist
    pd = fitdist(run_lengths', 'Exponential');
    lambda_hat = 1 / pd.mu;  % rate parameter λ = 1/mean

    % Optional: Plot histogram and fitted PDF
    % figure;
    % histogram(run_lengths, 'Normalization', 'pdf');
    % hold on;
    % x_vals = linspace(0, max(run_lengths)+1, 100);
    % y_vals = lambda_hat * exp(-lambda_hat * x_vals);
    % plot(x_vals, y_vals, 'r', 'LineWidth', 2);
    % xlabel('Run length');
    % ylabel('Probability Density');
    % title(sprintf('Exponential Fit (\\lambda = %.3f)', lambda_hat));
    % legend('Data', 'Exponential fit');
    % hold off;

function start_end_idx = find_run_regions_idx(binary_array, type, run_threshold)

    % Ensure input is a row vector and pad binary array with opposite
    % population
    binary_array = [mod(type+1,2) binary_array(:)' mod(type+1,2)];

    % Find transitions
    diff_array = diff(binary_array);
    if type
        start_idx = find(diff_array == 1);  % Start of runs
        end_idx = find(diff_array == -1); % End of runs
    else
        start_idx = find(diff_array == -1);  % Start of runs
        end_idx = find(diff_array == 1); % End of runs
    end

    start_end_idx(:,1) = start_idx'; start_end_idx(:,2) = end_idx';
    rem_idx = find( start_end_idx(:,2) - start_end_idx(:,1) < run_threshold);
    start_end_idx(rem_idx,:) = []; %remove too short runs
    

% function start_end_idx_by_percentages = sliding_window_regions(binary_array, type, window_sizes, sorted_thresholds)
%     % sliding_window_regions
%     % Input:
%     %   binary_array - Vector of 0s and 1s
%     %   type (bool) - 0 or 1 depending on if you want to know the runs of zeros or
%     %   ones
%     %   window_sizes - (1xN) array of desired stretches of window
%     %   sorted_thresholds - (1xM+1) normalized array of boundaries
%     % Output:
%     %   start_end_idx_by_percentages - (MxN) cell containing an (Lx2) array
%     %   with beginning and end
% 
%     % Ensure input is a row vector and pad binary array with opposite
%     % population
%     binary_array = [mod(type+1,2) binary_array(:)' mod(type+1,2)];
% 
%     % Find transitions
%     diff_array = diff(binary_array);
% 
% 
%     start_end_idx(:,1) = start_idx'; start_end_idx(:,2) = end_idx';
%     rem_idx = find( start_end_idx(:,2) - start_end_idx(:,1) < run_threshold);
%     start_end_idx(rem_idx,:) = []; %remove too short runs