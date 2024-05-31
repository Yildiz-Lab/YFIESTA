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
ptsin = []; ptsout = []; % this is so unnecessary, but for a ttest I have to do some stupid things
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
    [bin, bout, pin, pout, pts, num] = analyze_blips_vs(steptrace.data);
    blipin = [blipin; bin]; blipout = [blipout; bout];
    ptsin = [ptsin; pin]; ptsout = [ptsout; pout];
    totpts = totpts + pts; stepnum = stepnum + num;

    if mean(blipout) < -30
        fprintf(fname)
    end

    end

end

[mu,s1,s2] = beta_confidence(length(blipin),length(ptsin));
fprintf(strcat("Prob inside window (window size 3): ", num2str(round(mu,3)), " +/- [", num2str(round(s1,3)), ", ", num2str(round(s2,3)), "]", "\n"))
fprintf(strcat("Size ",num2str(round(-mean(blipin),2))," +/- ",num2str(round(std(blipin),2)),"\n"))

[mu0,s10,s20] = beta_confidence(length(blipout),length(ptsout));
fprintf(strcat("Prob outside window: ", num2str(round(mu0,3)), " +/- [", num2str(round(s10,3)), ", ", num2str(round(s20,3)), "]", "\n"))
fprintf(strcat("Size ",num2str(round(-mean(blipout),2))," +/- ",num2str(round(std(blipout),2)),"\n"))

fn1 = figure();
hold on
c_on = bar(1:2,[mu,mu0]);
er = errorbar(1:2,[mu,mu0],[s1-mu,s10-mu0],[s2-mu,s20-mu0],'Color',[0 0 0],'LineStyle','none','LineWidth',1);
ax = gca;
ax.YLim = [0,0.08];
ax.XLim = [0.3,2.7];
set(gcf,"Position",[250,250,200,300])

% ttest2(ptsin,ptsout)

% cdf
mdl_exp_cdf = fittype('real(gammainc(k*x,1))','indep','x');
xcdf = sort(dt(~isnan(dt)));
% xcdf = xcdf(xcdf>0.003);
ycdf = (1:length(xcdf))/length(xcdf);
cdffit = fit(xcdf',ycdf',mdl_exp_cdf,'start',[0.001])

f0 = figure();
histogram(blipin,'BinWidth',2,'Normalization','probability','DisplayName','By step')
hold on
histogram(blipout,'BinWidth',2,'Normalization','probability','DisplayName','Away from step')
legend('Location','northwest')
% set(gcf,"Position",[360,360,600,260])


f = figure();
subplot(1,2,1)
ax = gca;
hh = histogram(dt);
hh.BinWidth = 0.002;
ax.YLim = [0,90];
ax.XLim = [-0.002,0.024];
hold on
xlin = 0:0.0001:0.03;
plot(xlin,1.5*max(hh.Values)*exp(-cdffit.k*xlin),'LineWidth',2,'Color','r')
ax.LineWidth = 0.75; % Set the axes linewidth
ax.XColor = 'k'; % Set the color of x-axis
ax.YColor = 'k'; % Set the color of y-axis
set(ax, 'XColor', 'k', 'YColor', 'k'); % Set font size for axis labels and ticks, and set color of axes lines
set(findall(gca, 'Type', 'Text'), 'Color', 'k', 'FontName', 'Arial', 'FontSize', 10); % Set the color of all text objects to black

subplot(1,2,2)
hh = histogram(dx);
ax = gca;
hh.BinWidth = 2;
ax.YLim = [0,58];
ax.XLim = [-56,2];
ax.LineWidth = 0.75; % Set the axes linewidth
ax.XColor = 'k'; % Set the color of x-axis
ax.YColor = 'k'; % Set the color of y-axis
set(ax, 'FontName', 'Arial', 'FontSize', 10, 'XColor', 'k', 'YColor', 'k'); % Set font size for axis labels and ticks, and set color of axes lines
set(findall(gca, 'Type', 'Text'), 'Color', 'k', 'FontName', 'Arial', 'FontSize', 10); % Set the color of all text objects to black
set(gcf,"Position",[360,360,600,260])


% f2 = figure();
% % dt steps > 0
% subplot(2,2,1)
% ax = gca;
% hh = histogram(dt(step > 0));
% hh.BinWidth = 0.002;
% ax.YLim = [0,90];
% ax.XLim = [-0.002,0.024];
% 
% % dx steps > 0
% subplot(2,2,2)
% hh = histogram(dx(step > 0));
% ax = gca;
% hh.BinWidth = 2;
% ax.YLim = [0,58];
% ax.XLim = [-42,2];
% 
% % dt steps < 0
% subplot(2,2,3)
% ax = gca;
% hh = histogram(dt(step < 0));
% hh.BinWidth = 0.002;
% ax.YLim = [0,90];
% ax.XLim = [-0.002,0.024];
% 
% % dt steps < 0
% subplot(2,2,4)
% hh = histogram(dx(step < 0));
% ax = gca;
% hh.BinWidth = 2;
% ax.YLim = [0,58];
% ax.XLim = [-42,2];


f3 = figure();
% dt steps > 0
subplot(1,2,1)
hold on
ax = gca;
hh = histogram(dt(step > 0));
hh2 = histogram(dt(step < 0));
hh.BinWidth = 0.002; hh2.BinWidth = hh.BinWidth;
ax.YLim = [0,60];
ax.XLim = [-0.002,0.024];
ax.LineWidth = 0.75; % Set the axes linewidth
ax.XColor = 'k'; % Set the color of x-axis
ax.YColor = 'k'; % Set the color of y-axis
ax.LineWidth = 0.75; % Set the axes linewidth
ax.XColor = 'k'; % Set the color of x-axis
ax.YColor = 'k'; % Set the color of y-axis
set(ax, 'FontName', 'Arial', 'FontSize', 10, 'XColor', 'k', 'YColor', 'k'); % Set font size for axis labels and ticks, and set color of axes lines
set(findall(gca, 'Type', 'Text'), 'Color', 'k', 'FontName', 'Arial', 'FontSize', 10); % Set the color of all text objects to black

% dx steps > 0
subplot(1,2,2)
hold on
hh = histogram(dx(step > 0),'DisplayName','Before forward');
hh2 = histogram(dx(step < 0),'DisplayName','After backward');
ax = gca;
hh.BinWidth = 2; hh2.BinWidth = hh.BinWidth;
ax.YLim = [0,40];
ax.XLim = [-56,2];
legend('Location','Northwest','FontName', 'Arial', 'FontSize', 10)
ax.LineWidth = 0.75; % Set the axes linewidth
ax.XColor = 'k'; % Set the color of x-axis
ax.YColor = 'k'; % Set the color of y-axis
set(ax, 'FontName', 'Arial', 'FontSize', 10, 'XColor', 'k', 'YColor', 'k'); % Set font size for axis labels and ticks, and set color of axes lines
set(findall(gca, 'Type', 'Text'), 'Color', 'k', 'FontName', 'Arial', 'FontSize', 10); % Set the color of all text objects to black

set(gcf,"Position",[360,360,600,260])

end