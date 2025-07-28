function blip_percent = simple_blips_percent(fullfilename)

directory = [];

if isempty(directory)
    [filename, pathname] = uigetfile('*.mat','MultiSelect','on');
    % fullfilename = fullfile(pathname,filename);
    
    if ~isa(filename, 'cell')
        fullfilename{1} = fullfile(pathname,filename);
    else
        for i = 1:length(filename)
            fullfilename{i} = fullfile(pathname,filename{i});
        end
    end
else
    if ~isa(fullfilename, 'cell')
        fullfilename = {fullfilename};
    end
end

tot_steps = zeros(1,length(fullfilename));
blips = zeros(1,length(fullfilename));

for i = 1:length(fullfilename)

steptrace = load(fullfilename{i});

data = steptrace.data;

% JS edit because 2C data doesn't seem to have changepoints for some reason
% JS Edit 2025/03/07
% Get changepoints back in for 2C data with many Nans.
% Should have been careful changing the changepoint code.
nnidx = find(~isnan(data.trace(:,1)));

% on-axis
data.trace = data.trace(nnidx,:);
chp = find( abs(data.trace(2:end,3) - data.trace(1:end-1,3)) > 0);
data.trace(chp+1,5) = 1;

tot_steps(i) = sum(data.trace(:,5));
blips(i) = length(data.blips);

blip_percent = blips./tot_steps;
% blip_percent = sum(blips)/sum(tot_steps);


end