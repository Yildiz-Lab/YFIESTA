function [dt, dx, step] = analyze_blips_wrapper(directory)

%wrapper to convert x,y to polar coordinates and plot behavior for many
%different traces

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
    if contains(f(i).name,'._') %JS Edit to remove extra '._' that randomly show up sometimes
    f(i) = [];
    end
end
fnum = length(f);

dt = [];
dx = []; step = [];

for i=1:fnum
    
    if isfile(directory)
        [d,f,e] = fileparts(directory);
        % if contains(f,'._') %JS Edit to delete extra ._ from an error
        % f = f(3:end);
        % end
        steptrace = load(fullfile(d,strcat(f,e)));
    else
        fname = f(i).name;
        % if contains(fname,'._') %JS Edit to delete extra ._ from an error
        % fname = fname(3:end);
        % end
        steptrace = load(fullfile(cd,'/',fname));
    end
    
    if isfield(steptrace,'data')
    [dtprime, dxprime, stepprime] = analyze_blips(steptrace.data);
    dt = [dt, dtprime]; dx = [dx, dxprime]; step = [step, stepprime];
    end

end

f = figure();
subplot(1,2,1)
ax = gca;
hh = histogram(dt);
hh.BinWidth = 0.002;
ax.YLim = [0,90];
ax.XLim = [-0.002,0.024];


subplot(1,2,2)
hh = histogram(dx);
ax = gca;
hh.BinWidth = 2;
ax.YLim = [0,58];
ax.XLim = [-42,2];



end