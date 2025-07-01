%% Version 1 with just division, not the appropriate thing to do

% % Generate sample data (replace with actual data)
% data1 = randn(1, 1000);  % First dataset
% data2 = randn(1, 1000) + 1;  % Second dataset (shifted mean)
% 
% % Compute empirical CDFs
% [f1, x1] = ecdf(data1);
% [f2, x2] = ecdf(data2);
% x1 = x1(2:end); f1 = f1(2:end); x2 = x2(2:end); f2 = f2(2:end); %remove first because repeat
% 
% % Interpolate to get matching probabilities
% prob1 = interp1(x1, f1, sort(data1), 'linear', 'extrap');
% prob2 = interp1(x2, f2, sort(data1), 'linear', 'extrap'); % Match data1 to CDF of data2
% 
% % Normalize: Divide probability values
% normalizedValues = prob1 ./ prob2;
% normalizedValues = normalizedValues(normalizedValues >= min(data1));
% normalizedValues = normalizedValues(normalizedValues <= max(data1));
% 
% % Plot normalized histogram
% figure()
% numBins = 20;
% subplot(3,1,1)
% histogram(data1, numBins, 'FaceColor', 'b');
% subplot(3,1,2)
% histogram(data2, numBins, 'FaceColor', 'b');
% subplot(3,1,3)
% numBins = 20;
% histogram(normalizedValues, numBins, 'FaceColor', 'b');
% xlabel('Normalized Values');
% ylabel('Frequency');
% title('CDF-Normalized Data');


% % Generate sample data (replace with actual data)
% data1 = randn(1, 1000);  % First dataset (e.g., observed data)
% data2 = randn(1, 1000) + 0.5;  % Second dataset (e.g., prior with shifted mean)
% 
% % Define evaluation points (optional: use combined range of both datasets)
% evalPoints = linspace(min([data1, data2]), max([data1, data2]), 200);
% 
% % Estimate probability densities using kernel density estimation (KDE)
% pdf1 = ksdensity(data1, evalPoints);  % P_target(x)
% pdf2 = ksdensity(data2, evalPoints);  % P_prior(x)
% 
% % Compute weights: Normalize to avoid extreme values
% weights = pdf1 ./ pdf2;
% weights = weights / sum(weights);  % Normalize so sum of weights = 1
% 
% figure()
% plot(evalPoints, weights)
% 
% % Assign weights to data1 points using interpolation
% dataWeights = interp1(evalPoints, weights, data1, 'linear', 'extrap');
% 
% % Plot weighted histogram
% numBins = 20;
% histogram(data1, numBins, 'FaceColor', 'b', 'Normalization', 'pdf', ...
%           'DisplayStyle', 'bar', 'BinLimits', [min(data1), max(data1)], ...
%           'BinWidth', range(data1)/numBins, 'EdgeColor', 'k', 'FaceAlpha', 0.7);
% hold on;
% histogram(data1, numBins, 'Normalization', 'pdf', ...
%           'DisplayStyle', 'stairs', 'BinLimits', [min(data1), max(data1)], ...
%           'LineWidth', 2, 'EdgeColor', 'r');
% hold off;
% xlabel('Data1 Values');
% ylabel('Weighted Probability');
% title('Weighted Histogram of Data1');
% legend('Raw Data1 Histogram', 'Weighted Data1 Histogram');


% Generate sample data (replace with actual data)
data1 = randn(1, 1000);  % First dataset (e.g., observed data)
data2 = randn(1, 1000) + 0.75;  % Second dataset (e.g., prior with shifted mean)

data1 = hh.Data; %Actual data histogram, load "F:\MINFLUX JS\Dynein 2C\Compiled Traces\Stier\blips\behaving_blips\combined_test\head_separation_when_dips_occur_compiled.fig"
data2 = nn.Data; %Histogram of weighted head separation for normalization load "F:\MINFLUX JS\Dynein 2C\Compiled Traces\Stier\Stier_all_on-off-axis-sep.fig"

% Define bins for the histogram
numBins = 35;
binEdges = linspace(min([min(data1), min(data2)]), max([max(data1), max(data2)]), numBins+1);
binCenters = binEdges(1:end-1) + diff(binEdges)/2;  % Compute bin centers

% Compute probability densities at bin centers
pdf1 = ksdensity(data1, binCenters);  % P_target(x)
pdf2 = ksdensity(data2, binCenters);  % P_prior(x)

% Compute weights: Avoid division by zero
weights = pdf1 ./ pdf2;
weights(isnan(weights) | isinf(weights)) = 0;  % Replace NaN/Inf with zero

% Compute histogram counts for data1
[counts, ~] = histcounts(data1, binEdges);

% Apply weighting to histogram counts
weightedCounts = counts .* weights;

% Normalize to match probability density scale
weightedCounts = weightedCounts / sum(weightedCounts); 

% Plot histograms
figure;
bar(binCenters, counts / sum(counts), 'FaceColor', 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'k');  % Raw histogram
hold on;
bar(binCenters, weightedCounts, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'k');  % Weighted histogram
hold off;
xlabel('Data1 Values');
ylabel('Normalized Frequency');
title('Weighted Histogram');
legend('Raw Histogram', 'Weighted Histogram');

figure()
hold on
hh1 = histogram(data1, numBins, 'FaceColor', 'b', 'Normalization', 'pdf', ...
          'DisplayStyle', 'bar', 'BinLimits', [min(data1), max(data1)], ...
          'BinWidth', range(data1)/numBins, 'EdgeColor', 'k', 'FaceAlpha', 0.7);
hh2 = histogram(data2, numBins, 'FaceColor', 'r', 'Normalization', 'pdf', ...
          'DisplayStyle', 'bar', 'BinLimits', [min(data1), max(data1)], ...
          'BinWidth', range(data1)/numBins, 'EdgeColor', 'k', 'FaceAlpha', 0.5);

% hh1 = histogram(data1, numBins, 'FaceColor', 'b', ...
%           'DisplayStyle', 'bar', 'BinLimits', [min(data1), max(data1)], ...
%           'BinWidth', range(data1)/numBins, 'EdgeColor', 'k', 'FaceAlpha', 0.7);
% hh2 = histogram(data2, numBins, 'FaceColor', 'r', ...
%           'DisplayStyle', 'bar', 'BinLimits', [min(data1), max(data1)], ...
%           'BinWidth', range(data1)/numBins, 'EdgeColor', 'k', 'FaceAlpha', 0.5);

binCenters = 0.5*hh1.BinEdges(2:end) + 0.5*hh1.BinEdges(1:end-1);
normalized_values = hh1.Values(hh2.Values > 0) ./ hh2.Values(hh2.Values > 0);
plot(binCenters(hh2.Values > 0), normalized_values * max(hh1.Values) / max(normalized_values), 'k-', 'LineWidth', 2)
% first order error propagation means sigma_ratio = meanA / meanB * sqrt( (stdA / meanA)^2 + (stdB / meanB)^2 )
% since these are counting statistics, the ratio std/mu = 1/sqrt(N)
% This therefore simplifies to meanA/meanB * sqrt(1/meanA + 1/meanB)

xlabel("Head Separation at dips (nm)")
legend('Separation of heads at dip', 'Histogram of head separation', 'Normalized Probability of dip occuring');

figure()
uA = hh1.Values(hh2.Values > 0);
uB = hh2.Values(hh2.Values > 0);
sAB = uA ./ uB .* sqrt(1./uA + 1./uB);
errorbar(binCenters(hh2.Values > 0), normalized_values * max(hh1.Values) / max(normalized_values), sAB)