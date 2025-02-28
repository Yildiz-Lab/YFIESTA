% Make sure to load Molecule

mindistance = 300; %nm

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
hold on
scatter(distances, efo)
plot([min(distances);max(distances)], mean(efo)*ones(2,1),'r--')

figure()
hold on
scatter(distances, eco)
plot([min(distances);max(distances)], mean(eco)*ones(2,1),'r--')