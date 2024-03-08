function [m,c1,c2] = beta_confidence(a,b,conf)
% JS 2024/02/20

%Uses the beta distribution to calculate the error of percentage of a
%stochastic counting process with two options (i.e. stepping)
% https://youtu.be/8idr1WZ1A7Q?si=tYfVCQjD3X9CJ1IL
% https://youtu.be/ZA4JkHKZM50?si=1JOHh87X2fHG-Eoh
% https://youtu.be/juF3r12nM5A?si=w5s4e59eqjS4rTUc
% https://en.wikipedia.org/wiki/Beta_distribution

% Parameters:
    % a: Number of events of type a
    % b: Number of events not type a
    % conf: confidence interval desired
    
% Returns:
    % m: Most likely value of p in the binomial distribution for event a
    % (should be close to Na / (Na+Nb)
    % c1,c2: error based off of specified confidence interval (where cdf is
    % 0.5-conf/2 and 0.5+conf/2)

% make interpolation such that it is based off the total number of events
% as to make the correct resolution. Probably unnecessary will check both.

if nargin < 3
    conf = 0.95;
end

if conf > 1
    fprintf("Assuming percent, divide by 100 for confidence interval \n")
    conf = conf/100 %assume percent rather than fraction
end

x = 0:(a+b-1)^(-1):1;
% y = 0:0.001:1;

[~,idx1] = max(betapdf(x,a,b));
% [~,idx2] = max(betapdf(y,a,b));
m = (idx1-mod(length(x),2))/length(x);
% idx2 = (idx2-mod(length(y),2))/length(y);


[~,c1] = min(abs(betacdf(x,a,b)-(0.5-conf/2)));
[~,c2] = min(abs(betacdf(x,a,b)-(0.5+conf/2)));

c1 = c1/length(x);
c2 = c2/length(x);

end

