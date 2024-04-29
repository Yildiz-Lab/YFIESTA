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
blipin = []; blipout = [];
totpts = 0; stepnum = 0;

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
    [bin, bout, pts, num] = analyze_blips_vs(steptrace.data);
    blipin = [blipin; bin]; blipout = [blipout; bout];
    totpts = totpts + pts; stepnum = stepnum + num;

    if mean(blipout) < -30
        fprintf(fname)
    end

    end

end

fprintf(strcat("Prob inside window (window size 3): ", num2str(round(length(blipin)/(3*stepnum),2)),"\n"))
fprintf(strcat("Size ",num2str(round(-mean(blipin),2))," +/- ",num2str(round(std(blipin),2)),"\n"))
fprintf(strcat("Prob outside window: ", num2str(round(length(blipout)/(totpts-3*stepnum),2)),"\n"))
fprintf(strcat("Size ",num2str(round(-mean(blipout),2))," +/- ",num2str(round(std(blipout),2)),"\n"))

ttest2(blipin,blipout)

% cdf
mdl_exp_cdf = fittype('real(gammainc(k*x,1))','indep','x');
xcdf = sort(dt(~isnan(dt)));
ycdf = (1:length(dt(~isnan(dt))))/length(dt(~isnan(dt)));
cdffit = fit(xcdf',ycdf',mdl_exp_cdf,'start',[0.001])

f0 = figure();
histogram(blipin,'BinWidth',2,'Normalization','probability','DisplayName','By step')
hold on
histogram(blipout,'BinWidth',2,'Normalization','probability','DisplayName','Away from step')
legend('Location','northwest')

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