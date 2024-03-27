function plot_polar_conversion(r,theta,bindeg)

% Pass in radius (+,- no issue)
% Theta
% bindeg - number of degrees per bin

if nargin < 3
    bindeg = 15;
end

% plot in cool spherical histogram fashion
f = figure();
subplot(1,2,1)
polarscatter(theta,r)
set(gca,'ThetaZeroLocation','top')

subplot(1,2,2)
polarhistogram(theta,'BinWidth',bindeg*pi/180)
set(gca,'ThetaZeroLocation','top')

% subplot(2,2,3)
% polarhistogram(theta(r>0),'BinWidth',bindeg*pi/180)
% set(gca,'ThetaZeroLocation','top')
% 
% subplot(2,2,4)
% polarhistogram(theta(r<0),'BinWidth',bindeg*pi/180)
% set(gca,'ThetaZeroLocation','top')


end