function plot_polar_conversion(r,theta,binnum)

% Pass in radius (+,- no issue)
% Theta
% bindeg - number of degrees per bin

if nargin < 3
    binnum = 24;
end

% plot in cool spherical histogram fashion
f = figure();
subplot(1,2,1)
polarscatter(theta,r)
set(gca,'ThetaZeroLocation','top')


% theta_transform = theta_transform - pi*(theta_transform > pi);
fwd_mask = and(theta > -pi/2, theta < pi/2);

subplot(1,2,2)
polarhistogram(theta,'BinEdges',linspace(0,2*pi,binnum))
% polarhistogram(theta_transform,'BinEdges',linspace(0,2*pi,binnum))
hold on
ax = gca;
polarplot(median(theta(fwd_mask),'omitnan')*[1,1],[0,ax.RLim(2)],'r--')
polarplot(median(theta(~fwd_mask),'omitnan')*[1,1],[0,ax.RLim(2)],'r--')
fprintf(strcat("Mean fwd theta (deg): ", num2str(mean(theta(fwd_mask)*180/pi,'omitnan')), "\n"))
set(gca,'ThetaZeroLocation','top')

% subplot(2,2,3)
% polarhistogram(theta(r>0),'BinWidth',bindeg*pi/180)
% set(gca,'ThetaZeroLocation','top')
% 
% subplot(2,2,4)
% polarhistogram(theta(r<0),'BinWidth',bindeg*pi/180)
% set(gca,'ThetaZeroLocation','top')


end