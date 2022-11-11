function fStepStats()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global Config;

checkdir = Config.Directory{1};
actualdir = uigetdir(checkdir);

if length(Config.Time) > 1
    time = Config.Time(1);
else
    time = Config.Time;
end

StepandDwellHist_v2(actualdir, 0, time/1000);

end

