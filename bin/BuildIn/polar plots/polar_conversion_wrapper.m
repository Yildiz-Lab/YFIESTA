function [r,theta] = polar_conversion_wrapper(directory)

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

r = [];
theta = [];

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
    [rprime, thetaprime] = polar_conversion(steptrace.data,0);
    r = [r, rprime]; theta = [theta, thetaprime];
    end

end

plot_polar_conversion(r, theta, 24)

end