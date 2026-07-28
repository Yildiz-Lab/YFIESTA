function order = fStepSort(PathStats, ideal_mean_slope)
%fStepSort

% Description: sort paths by best possible traces to worst possible traces
%       to speed up analysis

% Parameters:
%       PathStats - PathStats structure as given in fPathStatsGui
%       ideal_mean_slope (optional) - the ideal slope for a stepping trace
%           to be partially used in ranking best to worst traces

% Returns:
%       Order - array of indices by ascending scores. The best score is the lowest.

if nargin < 2 %give a default reasonable stepping trace value
    ideal_mean_slope = 2.5;
end

nPath = size(PathStats,2);
score = zeros(1,nPath);

for n = 1:nPath
    frames = PathStats(n).Results(:,1); numdata = length(frames);
    onax = PathStats(n).PathData(:,4); offax = PathStats(n).PathData(:,5);
    
    % Sliding Window
    tau = [1, 5, 10, 20, 40, 80, 160, 320, 640];
    ntau = length(tau);
    tau = tau(tau < max(frames) - min(frames)); % make sure to only do if we can get a data point
    mu = zeros(1,length(tau)); sigma = zeros(1,length(tau));
    mu_confdiff = zeros(1,length(tau)); sigma_confdiff = zeros(1,length(tau));
    
    N = 0;
    % calculate each sliding window
    for j = 1:length(tau)
        dt = tau(j);
        
        vals = [];
        for i = 1:numdata
            [idx,~] = find(frames == frames(i) + dt);
            if ~isempty(idx)
                vals = [vals, onax(idx) - onax(i)];
            end
        end
        if isempty(vals)
            mu(j) = NaN; sigma(j) = NaN;
            mu_confdiff(j) = NaN; sigma_confdiff(j) = NaN;
        else
            % fit gaussian cdf to vals and get parameters
            mu(j) = mean(vals); sigma(j) = std(vals);
    
    %         mdl_norm_cdf = fittype('normcdf(x,mu,sigma)','indep','x');
    %         x = sort(vals);
    %         p = ((1:length(vals))-0.5)' ./ length(vals);
    %         
    %         % (x,p) should now be a cdf for a simple gaussian. Fit it
    %         ff = fit(x',p,mdl_norm_cdf,'start',[mu(j),sigma(j)]);
    %         mu(j) = ff.mu; sigma(j) = ff.sigma;
    %         
    %         cc = confint(ff, 0.95);
    %         mu_confdiff(j) = cc(2,1) - cc(1,1);
    %         sigma_confdiff(j) = cc(2,2) - cc(1,2);
        
            N = N + length(vals);
        end
        
    end
    
    % the better walkers will have the most linear increases in
    % mean and also increase (sqrt?) variance
    
    tau = tau(~isnan(mu));
    mu = mu(~isnan(mu)); sigma = sigma(~isnan(sigma));
    %mu_confdiff = mu_confdiff(~isnan(mu_confdiff)); sigma_confdiff = sigma_confdiff(~isnan(sigma_confdiff));
    
%     figure()
%     plot(tau, mu)
%     hold all
%     plot(tau, sigma)
%     plot(tau, mu_confdiff)
%     plot(tau, sigma_confdiff)
    
    % to create a score:
    %   mean(sigma)/sqrt(N*num_SW for this trace) + difference of slope from user expected / N^(inverse number of possible SW)
    % lowest score is the best score
    
    mpf = polyfit(tau, mu, 1); % do linear least squares fit on mu
    
    %score(n) = 1/sqrt(N)*(mpf(1) - ideal_mean_slope)^2;
    score(n) = N^(-1/ntau)*abs(mpf(1) - ideal_mean_slope) + (N*length(sigma))^(-1/2)*mean(sigma);
    %score(n) = N^(-1/ntau)*(abs(mpf(1) - ideal_mean_slope) + );
    %scoren(n) = numdata^(-1/ntau)*(mpf(1) - ideal_mean_slope)^2;
    
end

[scores, order] = sort(score);
%[scoresn, ordern] = sort(scoren)

end

