% quick sorting algorithm to quickly sort through traces based on the field
% of the matlab trace

% matlab traces have 

dc = uigetdir();

f = dir(fullfile(dc,'*.mat')); %JS Edit 220207
fnum = length(f);

for i=1:fnum

    fname = f(i).name;
    fname_cell{i} = f(i).name;
    [~,ONsteps{i},OFFsteps{i},dwells{i},dwells_for{i},dwells_back{i}] = StepandDwellHist_subfn(fullfile(dc,'/',fname),0,70);
    

end

% Now rank according to any criterion based on ONsteps, OFFsteps, dwells,
% dwells_for, dwells_back

A = OFFsteps;

[~,I] = sort(cellfun(@length,A));

% Now print off the names in the ascending order

for j=1:10
    fname_cell{I(j)}
end

%A = A(I);

