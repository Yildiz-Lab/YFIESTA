function drift_demonstration(Molecule)

% plot all molecules in their spatial coordinates with a time component by
% colormap

figure()
hold on

for j = 1:length(Molecule)
    t = Molecule(j).Results(:,2);
    t_norm = (t - min(t)) / (max(t) - min(t));
    scatter(Molecule(j).Results(:,3), Molecule(j).Results(:,4), 0.5, t_norm, 'filled')
    % scatter(Molecule(j).Results(:,3), Molecule(j).Results(:,4), 0.5, Molecule(j).Results(:,2)-Molecule(j).Results(1,2), 'filled')
    colormap("jet");
    % % Customize the color bar to show datetime instead of numeric values
    % num_ticks = 4;  % Number of ticks you want on the colorbar
    % tick_values = linspace(min(t_numeric), max(t_numeric), num_ticks);  % Get the tick positions in numeric form
    % % Convert the tick positions back to datetime
    % tick_labels = datestr(tick_values, 'yyyy-mm-dd HH:MM:SS');  % Customize the datetime format as needed
    % % Apply the tick labels to the colorbar
    % c.Ticks = tick_values;  % Set the tick positions on the colorbar
    % c.TickLabels = tick_labels; % Set the tick labels as datetime 
end


end