function [disp, sigma] = estimate_precision_error(directory)

% JS 2024/03/27

%wrapper to convert x,y to polar coordinates and plot behavior for many
%different traces. Check out polar_conversion.m for more details

%option to pass in directory directly or let user choose through gui

if nargin < 1
    directory = uigetdir()
end

if isfile(directory)
    fname = directory;
    fnum = 1;
else
    % Gather steps and dwells all in one folder
    cd = directory; %JS Edit 220207
    f = dir(fullfile(cd,'*.mat')); %JS Edit 220207
    fnum = length(f);
end

for i=fnum:-1:1
    if contains(f(i).name,'._') %JS Edit to ignore extra '._' that randomly show up sometimes MAC
    f(i) = [];
    end
end
fnum = length(f);

dr = [];
dt = [];

r = [];
t = [];

for i=1:fnum
    
    if isfile(directory)
        [d,f,e] = fileparts(directory);
        % if contains(f,'._') %JS Edit to delete extra ._ from an error
        % f = f(3:end);
        % end
        steptrace = load(fullfile(d,strcat(f,e)));
    else
        fname = f(i).name
        % if contains(fname,'._') %JS Edit to delete extra ._ from an error
        % fname = fname(3:end);
        % end
        steptrace = load(fullfile(cd,'/',fname));
    end
    
    if isfield(steptrace,'data')
        % convert to dist x^2 + y^2
        data = steptrace.data;
        dr = [dr; sqrt( (data.trace(1:end-1,1)-data.trace(2:end,1)).^2 + (data.trace(1:end-1,2)-data.trace(2:end,2)).^2 )];
        dt = [dt; data.time(2:end)-data.time(1:end-1)];
        r = [r; sqrt( (data.trace(1,1)-data.trace(end,1)).^2 + (data.trace(1,2)-data.trace(end,2)).^2 )];
        t = [t; data.time(end)-data.time(1)];
        % trace specific
        drind = sqrt(var(sqrt( (data.trace(1:end-1,1)-data.trace(2:end,1)).^2 + (data.trace(1:end-1,2)-data.trace(2:end,2)).^2 ),'omitnan'))/sqrt(2)
        dtind = mean(data.time(2:end)-data.time(1:end-1),'omitnan')
    end

end

sigma = sqrt(var(dr,'omitnan'))
sigma/sqrt(2)

meandt = mean(dt,'omitnan')

lifetime = mean(t, 'omitnan')

end