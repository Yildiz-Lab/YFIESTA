% function results = fit_dwell_cdfs(dwell)
% % fit_dwell_cdfs  Fit empirical CDF of dwell times with:
% %   1) Exponential CDF: F(t) = 1 - exp(-lambda*t)   (parameter lambda>0)
% %   2) Gamma (k=2, scale theta) CDF:
% %        F(t) = 1 - exp(-t/theta) * (1 + t/theta)  (theta>0)
% %
% % Usage:
% %   results = fit_dwell_cdfs(dwell)
% %
% % Returns struct with fields: x, ecdf, params_exp, params_gamma2,
% % gof_exp (R2,AIC,BIC,SSE), gof_gamma2, h (figure handle)
% 
% if nargin<1 || isempty(dwell)
%     error('Provide a nonempty vector of dwell times.');
% end
% 
% dwell = dwell(:);
% % Remove nonpositive values
% dwell = dwell(dwell>0);
% if isempty(dwell)
%     error('No positive dwell times.');
% end
% 
% % Empirical CDF using unique values (monotonic x)
% [x_ecdf, f_ecdf] = ecdf(dwell);
% % ecdf returns x starting from 0; remove duplicates with 0 probability if needed
% % Keep x>0 for fitting to avoid degenerate at zero
% mask = x_ecdf > 0;
% x = x_ecdf(mask);
% y = f_ecdf(mask);
% 
% n = numel(y); % number of points used in fit
% % plot(y,x)
% 
% % Model functions
% exp_cdf = @(lambda, t) 1 - exp(-lambda .* t);
% gamma2_cdf = @(theta, t) 1 - exp(-t./theta) .* (1 + t./theta);
% 
% % Objective SSE function (on y)
% sse_exp = @(logp) sum((y - exp_cdf(exp(logp), x)).^2);
% sse_gamma2 = @(logp) sum((y - gamma2_cdf(exp(logp), x)).^2);
% 
% % Initial guesses (on log scale)
% lambda0 = 1 / median(dwell);      % exp rate
% theta0  = median(dwell) / 1;      % scale guess
% 
% opts = optimset('Display','off','TolX',1e-8,'TolFun',1e-10,'MaxIter',1e6,'MaxFunEvals',1e6);
% 
% % Use fminsearch (does not require toolboxes). Fit log-parameters to enforce >0
% loglambda_hat = fminsearch(sse_exp, log(lambda0), opts);
% logtheta_hat  = fminsearch(sse_gamma2, log(theta0), opts);
% 
% lambda_hat = exp(loglambda_hat);
% theta_hat  = exp(logtheta_hat);
% 
% % Predicted CDFs at x
% y_exp = exp_cdf(lambda_hat, x);
% y_gamma2 = gamma2_cdf(theta_hat, x);
% 
% % Goodness-of-fit metrics
% res_exp = y - y_exp;
% res_gamma2 = y - y_gamma2;
% 
% SSE_exp = sum(res_exp.^2);
% SSE_gamma2 = sum(res_gamma2.^2);
% SST = sum((y - mean(y)).^2);
% 
% R2_exp = 1 - SSE_exp ./ SST;
% R2_gamma2 = 1 - SSE_gamma2 ./ SST;
% 
% % Number of fitted parameters
% k_exp = 1;
% k_gamma2 = 1;
% 
% % AIC and BIC (using residual variance SSE/n)
% AIC_exp = n * log(SSE_exp / n) + 2 * k_exp;
% AIC_gamma2 = n * log(SSE_gamma2 / n) + 2 * k_gamma2;
% 
% BIC_exp = n * log(SSE_exp / n) + k_exp * log(n);
% BIC_gamma2 = n * log(SSE_gamma2 / n) + k_gamma2 * log(n);
% 
% % Prepare results
% results.x = x;
% results.ecdf = y;
% results.params_exp = struct('lambda', lambda_hat);
% results.params_gamma2 = struct('theta', theta_hat);
% results.gof_exp = struct('R2', R2_exp, 'AIC', AIC_exp, 'BIC', BIC_exp, 'SSE', SSE_exp);
% results.gof_gamma2 = struct('R2', R2_gamma2, 'AIC', AIC_gamma2, 'BIC', BIC_gamma2, 'SSE', SSE_gamma2);
% 
% % Plot
% h = figure;
% hold on;
% % empirical CDF (step)
% stairs([0; y], [0; x], 'k', 'LineWidth', 1.5);
% % smooth plot domain
% tplot = linspace(0, max(y)*1.05, 200)';
% plot(tplot, exp_cdf(lambda_hat, tplot), 'r-', 'LineWidth', 1.4);
% plot(tplot, gamma2_cdf(theta_hat, tplot), 'b--', 'LineWidth', 1.4);
% legend('Empirical CDF', 'Exponential fit', 'Gamma (k=2) fit', 'Location','best');
% xlabel('Dwell time');
% ylabel('CDF');
% title('Empirical CDF and Fits');
% grid on;
% hold off;
% results.h = h;
% 
% % Display summary
% fprintf('Exponential fit: lambda = %.6g\n', lambda_hat);
% fprintf('  R^2 = %.6f, AIC = %.4f, BIC = %.4f, SSE = %.6g\n', R2_exp, AIC_exp, BIC_exp, SSE_exp);
% fprintf('Gamma (k=2) fit: theta = %.6g\n', theta_hat);
% fprintf('  R^2 = %.6f, AIC = %.4f, BIC = %.4f, SSE = %.6g\n', R2_gamma2, AIC_gamma2, BIC_gamma2, SSE_gamma2);
% 
% end

function results = fit_dwell_cdfs(dwell)
% fit_dwell_cdfs  Fit dwell times with:
%   - Exponential (MLE via fitdist)
%   - Gamma with fixed shape k=2 (MLE for scale theta)
% Returns struct with params, gof (R2,AIC,BIC,logL), and figure handle.

if nargin<1 || isempty(dwell)
    error('Provide a nonempty vector of dwell times.');
end

dwell = dwell(:);
nData = numel(dwell);
if nData==0, error('No positive dwell times.'); end

% Empirical CDF (unique sorted values, final value = 1)
xs = sort(dwell);
[xu, ~, ic] = unique(xs);
counts = accumarray(ic,1);
ecdf_y = cumsum(counts) / nData;
x = xu;
y = ecdf_y;  % y(end) == 1

% Fit exponential using fitdist (MLE). fitdist returns mean 'mu'.
exp_pd = fitdist(dwell,'Exponential');
mu_hat = exp_pd.mu;
lambda_hat = 1 / mu_hat;

% Exponential log-likelihood
logL_exp = sum(log(lambda_hat) - lambda_hat .* dwell);

% Fit gamma(k=2) by MLE for theta (scale). PDF for k=2: f(t) = t * exp(-t/theta) / theta^2
negloglik_gamma2 = @(logtheta) -sum( log(dwell) - dwell./exp(logtheta) - 2*logtheta );
% initial guess: theta0 = mean/2
logtheta0 = log(mean(dwell)/2);
opts = optimset('TolX',1e-12,'MaxIter',1e4,'Display','off');
logtheta_hat = fminsearch(negloglik_gamma2, logtheta0, opts);
theta_hat = exp(logtheta_hat);
logL_gamma2 = -negloglik_gamma2(logtheta_hat);

% Predicted CDFs
exp_cdf = @(t) 1 - exp(-lambda_hat .* t);
gamma2_cdf = @(t) 1 - exp(-t./theta_hat) .* (1 + t./theta_hat);  % analytic for k=2

% Evaluate fits on CDF x points
y_exp = exp_cdf(x);
y_gamma2 = gamma2_cdf(x);

% Goodness-of-fit: R^2 on CDF points (using SST of empirical CDF)
SSE_exp = sum((y - y_exp).^2);
SSE_gamma2 = sum((y - y_gamma2).^2);
SST = sum((y - mean(y)).^2);
R2_exp = 1 - SSE_exp./SST;
R2_gamma2 = 1 - SSE_gamma2./SST;

% AIC and BIC using log-likelihood
k_exp = 1;        % one fitted parameter (lambda)
k_gamma2 = 1;     % one fitted parameter (theta)

AIC_exp = 2*k_exp - 2*logL_exp;
AIC_gamma2 = 2*k_gamma2 - 2*logL_gamma2;

BIC_exp = k_exp*log(nData) - 2*logL_exp;
BIC_gamma2 = k_gamma2*log(nData) - 2*logL_gamma2;

% Output struct
results.x = x;
results.ecdf = y;
results.params_exp = struct('mu', mu_hat, 'lambda', lambda_hat);
results.params_gamma2 = struct('theta', theta_hat, 'k', 2);
results.gof_exp = struct('R2', R2_exp, 'AIC', AIC_exp, 'BIC', BIC_exp, 'logL', logL_exp, 'SSE', SSE_exp);
results.gof_gamma2 = struct('R2', R2_gamma2, 'AIC', AIC_gamma2, 'BIC', BIC_gamma2, 'logL', logL_gamma2, 'SSE', SSE_gamma2);

% Plot empirical CDF and fits
h = figure;
hold on;
% step empirical CDF
stairs([0; x], [0; y], 'w', 'LineWidth', 1.5);
tplot = linspace(0, max(x)*1.05, 300)';
plot(tplot, exp_cdf(tplot), 'r-', 'LineWidth', 1.4);
plot(tplot, gamma2_cdf(tplot), 'b--', 'LineWidth', 1.4);
legend('Empirical CDF', 'Exponential fit', 'Gamma (k=2) fit', 'Location','best');
xlabel('Dwell time');
ylabel('CDF');
title('Empirical CDF and Fits');
grid on;
hold off;
results.h = h;

% Display summary
fprintf('Exponential (MLE): lambda = %.6g (mu = %.6g)\n', lambda_hat, mu_hat);
fprintf('  R^2 = %.6f, logL = %.4f, AIC = %.4f, BIC = %.4f\n', R2_exp, logL_exp, AIC_exp, BIC_exp);
fprintf('Gamma (k=2) (MLE): theta = %.6g\n', theta_hat);
fprintf('  R^2 = %.6f, logL = %.4f, AIC = %.4f, BIC = %.4f\n', R2_gamma2, logL_gamma2, AIC_gamma2, BIC_gamma2);

end
