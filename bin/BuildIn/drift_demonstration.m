function drift_demonstration(Molecule)

% plot all molecules in their spatial coordinates with a time component by
% colormap

figure()
hold on

for j = 1:length(Molecule)
    t = Molecule(j).Results(:,2);
    if length(t) > 300
        t_norm = (t - min(t)) / (max(t) - min(t));
        scatter(Molecule(j).Results(:,3), Molecule(j).Results(:,4), 0.5, t_norm, 'filled')
        % scatter(Molecule(j).Results(:,3), Molecule(j).Results(:,4), 0.5, Molecule(j).Results(:,2)-Molecule(j).Results(1,2), 'filled')
        colormap("jet");
    end
end

colorbar
set(gca, 'YDir','reverse')