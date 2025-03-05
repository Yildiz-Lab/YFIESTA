function Plot2CStepStats(xydiff, xy_deltatxy_step)
%Parameters
% xydiff (Nx2 array) - contains all data on the separation of both channels in x and y
% xy_deltatxy_step (nx7 array) - each containing information about a changepoint step with the following format:
%   [channel, normalized trace time (s), x_head_separation (nm), y_head_separation (nm), deltat (s), deltax (nm), deltay (nm)];

% Returns:
% 1. Histograms of xy-trace difference in two color separation
% 2. Scatter plots of step size and dwell time depending on interhead separation

%% Individual plot options

figure()
% On-axis Inter-head separation
subplot(1,2,1)
histogram(xydiff(:,1))
title('On-axis Separation (nm)')

% Off-axis Inter-head separation
subplot(1,2,2)
histogram(xydiff(:,2))
title('Off-axis Separation (nm)')

figure()
% Step size
subplot(3,1,1)
scatter(xy_deltatxy_step(:,3), xy_deltatxy_step(:,6), 50, 'k', 'filled')
ylabel('\Delta x')
title('Dependence on inter-head separation')

subplot(3,1,2)
scatter(xy_deltatxy_step(:,4), xy_deltatxy_step(:,7), 50, 'k', 'filled')
ylabel('\Delta x')
title('Dependence on inter-head separation')

% Dwell
subplot(3,1,3)
scatter(xy_deltatxy_step(:,3), xy_deltatxy_step(:,5), 50, 'k', 'filled')
ylabel('Dwell time')
xlabel('Head Separation (nm)')

end