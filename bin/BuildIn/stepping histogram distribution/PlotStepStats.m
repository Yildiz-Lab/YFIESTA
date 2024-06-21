function PlotStepStats(tracenum, onsteps, offsteps, dwells, dwells_for, dwells_back, options, savename)
%Plot the characteristic statistics in nice figures

if nargin < 7
    options.Poissonk1 = 1;
    options.Poissonk2 = 1;
    options.DoubleExp = 1;
    options.FwdDwells = 1;
    savename = [];
elseif nargin < 8
    savename = [];
end

%JS Edit 220213
% Some general statistics that I would like to show in the plot legend
% because I know Ahmet will ask about it
N = sum(~isnan(onsteps)); %Ignore nans
Nfor = length(onsteps(onsteps > 0));
Nback = length(onsteps(onsteps < 0));

h =  findobj('type','figure');
n = length(h);

% On-axis step histogram
figure(n+1)
subplot(2,3,1)
hh = histogram(onsteps,'BinWidth',1.5);
axis([-40,48,0,40]);
set(gca, 'XTick', [-80 -32 -24 -16 -8 0 8 16 24 32 40 100]);
set(gca, 'ylim', [0 max(hh.Values)+5])
xlabel('step size (nm)');
title ('on-axis steps')
legend(sprintf(' traces = %.0f \n N = %.0f \n forward = %.0f \n backward = %.0f', [tracenum, N, Nfor(1), Nback(1)]))

% use betacdf to get the 95% confidence intervals
totalsteps = length(onsteps)+length(offsteps);
[m,c1,c2] = beta_confidence(Nback,totalsteps-Nback);
fprintf(strcat("Backwards / total stepping ", num2str(round(m,4)), " (", num2str(round(c1,4)), ",", num2str(round(c2,4)), ") \n"))
[m,c1,c2] = beta_confidence(Nback,Nfor);
fprintf(strcat("Backwards / forward stepping ", num2str(round(m,4)), " (", num2str(round(c1,4)), ",", num2str(round(c2,4)), ") \n"))

% Off-axis step histogram
subplot(2,3,2)
histogram(offsteps,'BinWidth',1.5);
axis([-32,32,0,50]);
set(gca, 'XTick', [-30 -25 -20 -15 -10 -5 0 5 10 15 20 25 30]);
% ylim([0 70])
xlabel('step size (nm)');
title ('off-axis steps')

mdl_gauss = fittype('normcdf(x,mu,sigma)','indep','x');
% mdl_gauss2 = fittype('A*normcdf(x,mu1,7.)+(1-A)*normcdf(x,mu2,7.)','indep','x');
X = sort(abs(offsteps));
Y = linspace(0,1,length(X));
try
fittedmdl = fit(X,Y',mdl_gauss,'start',[8.,8.])
legend(sprintf(' N = %.0f \n mean = %.3f \n std = %.3f', [length(offsteps), fittedmdl.mu, fittedmdl.sigma]))
catch
    legend(sprintf(' N = %.0f', [length(offsteps)]))
end
% fittedmdl2 = fit(X,Y',mdl_gauss2,'start',[0.7,8.,20.])



% use betacdf to get the 95% confidence intervals
[m,c1,c2] = beta_confidence(length(offsteps),totalsteps-length(offsteps));
fprintf(strcat("Side / total stepping ", num2str(round(m,4)), " (", num2str(round(c1,4)), ",", num2str(round(c2,4)), ") \n"))
[m,c1,c2] = beta_confidence(length(offsteps),Nfor);
fprintf(strcat("Side / forward stepping ", num2str(round(m,4)), " (", num2str(round(c1,4)), ",", num2str(round(c2,4)), ") \n"))

% Dwell histogram
% JS 220207 IDK who thought setting the xlimits to 0.25 seconds was a good
% idea, but getting rid of that allows us to see the plots.
subplot(2,3,3)
% histogram(dwells);
histogram(dwells, 'BinWidth', 0.01);
set(gca,'XLim',[-0.05,0.8])
xlabel('Time (s)');
ylabel('Counts');
title ('dwell times')

% Forward_Dwell histogram
subplot(2,3,4)
% histogram(dwells_for);
histogram(dwells_for,'BinWidth',0.01);
set(gca,'XLim',[-0.05,0.8])
xlabel('Time (s)');
title ('forward dwell times')

% BackwardDwell histogram
subplot(2,3,6)
% histogram(dwells_back);
histogram(dwells_back,'BinWidth',0.01);
set(gca,'XLim',[-0.05,0.8])
xlabel('Time (s)');
title ('backward dwell times')

%% JS Edit 220307
% make a CDF plot and fit
subplot(2,3,5)
if options.FwdDwells
    obj0 = cdfplot(dwells_for);
else
    obj0 = cdfplot(dwells);
end
xlabel('Time (s)');
hold on

% Define model according to https://www.mathworks.com/matlabcentral/answers/850245-how-to-do-curve-fitting-by-a-user-defined-function
% mdl_gamma_pdf = fittype('A * x*k^2*exp(-k*x)','indep','x');
mdl_exp_cdf = fittype('real(gammainc(k*x,1))','indep','x'); %single exponential
mdl_gamma_cdf = fittype('real(gammainc(k*x,2))','indep','x'); %k=2 fixed
% also should work: fittype(gamcdf(x,2,1/k))
% mdl_gamma_mixed_cdf = fittype('k1*k2/(k1+k2)*(exp((-k1*x)+exp(-k2*x)))','indep','x');
mdl_de_cdf = fittype('1-((p)*exp(-a*x) + (1-p)*exppdf(-b*x))','indep','x');
X = obj0.XData(2:end-1); % get rid of infs
Y = obj0.YData(2:end-1);
% fittedmdlmix = fit(X',Y',mdl_gamma_mixed_cdf,'start',[3.])
xsortminusnan = sort(dwells);
xsortminusnan = xsortminusnan(~isnan(xsortminusnan));
n = size(xsortminusnan,1);
p = ((1:n)-0.5)' ./ n;
ylog = -log(1 - p);
muHat = ylog \ xsortminusnan
muMLE = expfit(xsortminusnan)

if options.Poissonk1
    [fittedmdl,gof] = fit(X',Y',mdl_exp_cdf,'start',[3.])
    cfint = confint(fittedmdl);
plot(X',fittedmdl(X'),'r--','DisplayName',"Fitted 1/theta = "+num2str(round(fittedmdl.k,1))+" ("+num2str(round(cfint(1),1))+","+num2str(round(cfint(2),1))+") R2="+num2str(round(gof.rsquare,2)))
end
if options.Poissonk2
    [fittedmdl2,gof2] = fit(X',Y',mdl_gamma_cdf,'start',[3.])
    cfint2 = confint(fittedmdl2);
plot(X',fittedmdl2(X'),'k--','DisplayName',"Fitted 1/theta2 = "+num2str(round(fittedmdl2.k,1))+" ("+num2str(round(cfint2(1),1))+","+num2str(round(cfint2(2),1))+") R2="+num2str(round(gof2.rsquare,2)))
end
if options.DoubleExp
    % double exponential is 1 - (p*exp(-a*x)+(1-p)*exp(-b*x))
    % de_pdf = @(x,p,a,b) (p)*exppdf(x,a) + (1-p)*exppdf(x,b);
    opts = fitoptions('Method', 'NonlinearLeastSquares', ...
                     'Lower',[0. 0. 0.], 'Upper',[200. 200. 1.]);
    [fittedmdl3,gof3] = fit(X',Y',mdl_de_cdf,opts) %,'start',[3.,3.,0.5])
    % mle(dwells,'pdf',de_pdf,'start',[.5,-2,-3])
    plot(X',fittedmdl3(X'),'b--','DisplayName',"Fitted p = "+num2str(fittedmdl3.p))
end
mdl_glue_cdf = fittype('1-(p*exp(-a*x) + (1-p)*real(gammainc(b*x,2)))','indep','x')
opts = fitoptions('Method', 'NonlinearLeastSquares', ...
                     'Lower',[0. 15. 0.], 'Upper',[200. 30. 1.]);
[fittedmdl4,gof4] = fit(X',Y',mdl_glue_cdf,opts) %,'start',[3.,3.,0.5])

legend()
% set(gca,'XLim',[0,2])

if options.FwdDwells
    subplot(2,3,4)
else
    subplot(2,3,3)
end
hold on
ax = gca;
children = ax.Children;

tt = linspace(0,max(children.BinEdges),500);

% mdl_cdf = fittype('A*exp(-k*x)','indep','x');
if options.Poissonk1
    k = fittedmdl.k;
plot(tt, max(children.Values)*exp(-k.*tt), 'r--'); %single exponential
end
if options.Poissonk2
    k = fittedmdl2.k;
plot(tt, max(children.Values)/max(k^2/2.*tt.*exp(-k.*tt))*k^2/2.*tt.*exp(-k.*tt), 'k--'); %double exponential
end

if ~isempty(savename)
    savefig(savename)
end

end

