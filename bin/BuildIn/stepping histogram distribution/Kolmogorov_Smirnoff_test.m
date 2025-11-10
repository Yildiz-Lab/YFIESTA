% Kolmogorov-Smirnoff Test

% input data1 and data2 and use it to determine whether the samples are
% different or the same distribution

% Run combinations from workspace (Look for
% KSTest_stepping_hist_D1848E_velocity_pvals.mat)

datastruct = {data08, data20, data50, data1000};
pvals = zeros(length(datastruct), length(datastruct));

% to avoid number problems, we should bootstrap just a subset and see how
% it looks

num_elements = 400;
% num_elements = max([length(data08), length(data20), length(data50), length(data1000)]);

% p-tests lie if N of samples are not close to each other. Therefore, we
% randomly bootstrap 
% And let's do multiple iterations and take the average

for k = 1:10000
    randidx = zeros(num_elements, length(datastruct));
    for i = 1:length(datastruct)
        randidx(:,i) = randi(length(datastruct{i}), 1, num_elements);
    end

    for i = 1:length(datastruct)
        for j = 1:length(datastruct)
            if j - i > 0
                data1 = datastruct{i};
                data2 = datastruct{j};
                [h,p] = kstest2(data1(randidx(:,i)), data2(randidx(:,j)));
                % calculate a running average
                pvals(i,j) = ((k-1)*pvals(i,j) + p)/k;
            end
        end
    end

end

pvals

pvals_utest = zeros(length(datastruct), length(datastruct));

for k = 1:10000
    randidx = zeros(num_elements, length(datastruct));
    for i = 1:length(datastruct)
        randidx(:,i) = randi(length(datastruct{i}), 1, num_elements);
    end

    for i = 1:length(datastruct)
        for j = 1:length(datastruct)
            if j - i > 0
                data1 = datastruct{i};
                data2 = datastruct{j};
                [p,h] = ranksum(data1(randidx(:,i)), data2(randidx(:,j)));
                % [p,h] = ranksum(data1, data2);
                % calculate a running average
                pvals_utest(i,j) = ((k-1)*pvals_utest(i,j) + p)/k;
            end
        end
    end

end

pvals_utest



% Also Gmm fit
X = vertcat(data08, data20, data50, data1000);
% X = data20;

options = statset('MaxIter',1000);
gm = fitgmdist(X, 8, 'Options', options);

figure()
histogram(X,'BinWidth',2,'Normalization','pdf');
hold on
for i = 1:length(gm.mu)
    hold on
    x = -60:0.1:80;
    % y = 1.5*gm.ComponentProportion(i) * normpdf(x,gm.mu(i),gm.Sigma(1,1,i));
    y = gm.ComponentProportion(i) * normpdf(x,gm.mu(i),4);
    plot(x, y);

end

for i = 1:length(gm.mu)
    hold on
% Create a unique label for each point
    % label_text = sprintf(num2str(round(gm.mu,1)));
    label_text = num2str(round(gm.mu(i),1));

    % Add the text to the plot
    % text(gm.mu(i), gm.ComponentProportion(i), label_text, 'FontSize', 8, 'Color', 'red', ...
    %      'HorizontalAlignment', 'center'); 
    text(gm.mu(i), 0.04, label_text, 'FontSize', 10, 'Color', 'red', ...
         'HorizontalAlignment', 'center');
end
% text(gm.mu, gm.ComponentProportion, centers)
% Single overlay
gm1 = fitgmdist(X, 1, 'Options', options);
x = -60:0.1:80;
% y = 1.5*gm.ComponentProportion(i) * normpdf(x,gm.mu(i),gm.Sigma(1,1,i));
y = normpdf(x,mean(X),sqrt(var(X)));
plot(x, y);



% ttest2

datastruct = {WT, PO4, D1848E};
pvalsttest = zeros(length(datastruct), length(datastruct));

% to avoid number problems, we should bootstrap just a subset and see how
% it looks

num_elements = 800;

% p-tests lie if N of samples are not close to each other. Therefore, we
% randomly bootstrap 
% And let's do multiple iterations and take the average

for k = 1:10000
    randidx = zeros(num_elements, length(datastruct));
    for i = 1:length(datastruct)
        randidx(:,i) = randi(length(datastruct{i}), 1, num_elements);
    end

    for i = 1:length(datastruct)
        for j = 1:length(datastruct)
            if j - i > 0
                data1 = datastruct{i};
                data2 = datastruct{j};
                [h,p] = ttest2(data1(randidx(:,i)), data2(randidx(:,j)));
                % calculate a running average
                pvalsttest(i,j) = ((k-1)*pvalsttest(i,j) + p)/k;
            end
        end
    end

end

pvalsttest