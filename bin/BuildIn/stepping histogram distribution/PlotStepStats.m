function PlotStepStats(onsteps, offsteps, dwells, dwells_for, tracenum)
%Plot the characteristic statistics in nice figures

%JS Edit 220213
% Some general statistics that I would like to show in the plot legend
% because I know Ahmet will ask about it
N = length(onsteps);
Nfor = length(onsteps(onsteps > 0));
Nback = length(onsteps(onsteps < 0));

h =  findobj('type','figure');
n = length(h);

% On-axis step histogram
figure(n+1)
histogram(onsteps,'BinWidth',1);
axis([-40,48,0,40]);
set(gca, 'XTick', [-80 -32 -24 -16 -8 0 8 16 24 32 40 100]);
ylim([0 70])
xlabel('On-axis step size (nm)');
ylabel('Counts');
title ('on-axis step size distribution')
legend(sprintf(' traces = %.0f \n N = %.0f \n forward = %.0f \n backward = %.0f', [tracenum, N, Nfor(1), Nback(1)]))

% Off-axis step histogram
figure(n+2)
histogram(offsteps,'BinWidth',1);
axis([-40,48,0,100]);
set(gca, 'XTick', [-80 -32 -24 -16 -8 0 8 16 24 32 40 100]);
ylim([0 70])
xlabel('Off-axis step size (nm)');
ylabel('Counts');
title ('off-axis step size distribution')

% Dwell histogram
% JS 220207 IDK who thought setting the xlimits to 0.25 seconds was a good
% idea, but getting rid of that allows us to see the plots.
figure(n+3)
%histogram(dwells,'BinWidth',.001);
histogram(dwells, 'BinWidth', 1.5);
%xlim([0 0.25])
%ylim([0 50])
xlabel('Time (sec)');
ylabel('Counts');
title ('dwell time distribution')
% Forward_Dwell histogram
figure(n+4)
histogram(dwells_for,'BinWidth',1.5);
%histogram(dwells, 'BinWidth', 0.3);
%xlim([0 0.25])
%ylim([0 30])
xlabel('Time (sec)');
ylabel('Counts');
title ('forward dwell time distribution')
% BackwardDwell histogram
figure(n+5)
% histogram(dwells_back,'BinWidth',.001);
% %histogram(dwells, 'BinWidth', 0.3);
% %xlim([0 0.25])
% ylim([0 20])
% xlabel('Time (sec)');
% ylabel('Counts');
% title ('backward dwell time distribution')

%% JS Edit 220307
% make a CDF plot and fit
obj0 = cdfplot(dwells_for);
hold on

% Define model according to https://www.mathworks.com/matlabcentral/answers/850245-how-to-do-curve-fitting-by-a-user-defined-function
% mdl_gamma_pdf = fittype('A * x*k^2*exp(-k*x)','indep','x');
mdl_gamma_cdf = fittype('real(gammainc(k*x,2))','indep','x');
X = obj0.XData(2:end-1); % get rid of infs
Y = obj0.YData(2:end-1);
fittedmdl = fit(X',Y',mdl_gamma_cdf,'start',[3.])

plot(fittedmdl)
legend("Fitted k = "+num2str(fittedmdl.k))

figure(n+4)
hold on
ax = gca;
children = ax.Children;
A = length(children.Data);

k = fittedmdl.k;
tt = linspace(0,max(children.BinEdges),100);
plot(tt, A*k^2.*tt.*exp(-k.*tt));

% tt
% k^2.*tt.*exp(-k.*tt)


end

