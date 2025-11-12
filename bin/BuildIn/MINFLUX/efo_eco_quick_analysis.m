% Make sure to load Molecule

mindistance = 100; %nm

distances = [];
efo = [];
eco = [];

for i = 1:length(Molecule)
    
    distances = [distances; Molecule(i).Results(end,6)];
    efo = [efo; Molecule(i).Results(end,7)];
    eco = [eco; Molecule(i).Results(end,8)];

end

distances = distances(distances > mindistance);
efo = efo(distances > mindistance);
eco = eco(distances > mindistance);

figure()
subplot(1,2,1)
hold on
scatter(distances, efo)
plot([min(distances);max(distances)], mean(efo)*ones(2,1),'r--')
subplot(1,2,2)
histogram(efo)

figure()
subplot(1,2,1)
hold on
scatter(distances, eco)
plot([min(distances);max(distances)], mean(eco)*ones(2,1),'r--')
subplot(1,2,2)
histogram(efo)


%% precision error
dr = [];
dt = [];

r = [];
t = [];

mean_sigma = [];

for i = 1:length(Molecule)
    
    if Molecule(i).Results(end,6) > mindistance
        dr = [dr; sqrt( (Molecule(i).Results(1:end-1,3)-Molecule(i).Results(2:end,3)).^2 + (Molecule(i).Results(1:end-1,2)-Molecule(i).Results(2:end,2)).^2 )];
        dt = [dt; Molecule(i).Results(2:end,2)-Molecule(i).Results(1:end-1,2)];
        % r = [r; sqrt( (data.trace(1,1)-data.trace(end,1)).^2 + (data.trace(1,2)-data.trace(end,2)).^2 )];
        % t = [t; data.time(end)-data.time(1)];
        % % trace specific
        % drind = sqrt(var(sqrt( (data.trace(1:end-1,1)-data.trace(2:end,1)).^2 + (data.trace(1:end-1,2)-data.trace(2:end,2)).^2 ),'omitnan'))/sqrt(2)
        % dtind = mean(data.time(2:end)-data.time(1:end-1),'omitnan')

        % Let's also calculate just raw standard error based on trace fits.
        % (Remember that the trace fits need to be accurate)

        % mean_sigma = [mean_sigma, std(sqrt((data.trace(:,1)-data.trace(:,3)).^2 + (data.trace_yx(:,1)-data.trace_yx(:,3)).^2),'omitnan')];
    end

end

sigma = sqrt(var(dr,'omitnan'))
sigma/sqrt(2)

meandt = mean(dt,'omitnan')