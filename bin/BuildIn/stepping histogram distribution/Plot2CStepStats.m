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
x = xy_deltatxy_step(:,3);
y = xy_deltatxy_step(:,6);
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
histogram(x_on(leading_idx,2),'BinWidth',2)
title('Leading head')
ax=gca; set(ax,'XLim',long_axis_range);
set(ax, 'FontName', 'Arial', 'FontSize', 10, 'TickDir', 'out', 'LineWidth', 1, 'Box', 'off', 'XColor', 'k', 'YColor', 'k');

subplot(2,2,ceil(1+flip_reference_head/2))
histogram(x_on(trailing_idx,2),'BinWidth',2)
title('Trailing head')
ax=gca; set(ax,'XLim',long_axis_range);
set(ax, 'FontName', 'Arial', 'FontSize', 10, 'TickDir', 'out', 'LineWidth', 1, 'Box', 'off', 'XColor', 'k', 'YColor', 'k');

subplot(2,2,ceil(3-flip_reference_head/2))
histogram(x_off(right_idx,2),'BinWidth',2)
title('Right head')
ax=gca; set(ax,'XLim',short_axis_range);
set(ax, 'FontName', 'Arial', 'FontSize', 10, 'TickDir', 'out', 'LineWidth', 1, 'Box', 'off', 'XColor', 'k', 'YColor', 'k');

subplot(2,2,ceil(3+flip_reference_head/2))
histogram(x_off(left_idx,2),'BinWidth',2)
title('Left head')
ax=gca; set(ax,'XLim',short_axis_range);
set(ax, 'FontName', 'Arial', 'FontSize', 10, 'TickDir', 'out', 'LineWidth', 1, 'Box', 'off', 'XColor', 'k', 'YColor', 'k');


% And now the fraction of steps taken by each population, at least in
% on-axis

figure('Position',[311 266 722 420])
subplot(1,2,1)
% binedges = -95:10:95; %make symmetric across 0 for ease of flipping
binedges = -60:8:60; %make symmetric across 0 for ease of flipping
last_bin_show = 5;

hold on
h_trail = histogram(x_on(trailing_idx,2),'BinEdges',binedges,'DisplayName','Trailing');
h_lead = histogram(x_on(leading_idx,2),'BinEdges',binedges,'DisplayName','Leading');
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

% 
% % --- Helper function to keep parameters in bounds ---
% function nll = penalizedNegLL(p, nllFunc)
%     if p(1) < 0 || p(1) > 1 || p(2) <= 0 || p(3) <= 0
%         nll = Inf;
%     else
%         nll = nllFunc(p);
%     end
