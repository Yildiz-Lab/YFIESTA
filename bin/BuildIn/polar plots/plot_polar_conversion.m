function plot_polar_conversion(r,theta,binnum)

% JS 2024/03/27

% Description
% plot polar version of steps, the histogram, and a heat map

% Parameters
% Pass in radius (+,- no issue)
% Theta
% bindeg - number of degrees per bin

% 

if nargin < 3
    binnum = 24;
end

% plot in cool spherical histogram fashion
f = figure();
subplot(1,3,1)
polarscatter(theta,r)
set(gca,'ThetaZeroLocation','top')


% theta_transform = theta_transform - pi*(theta_transform > pi);
fwd_mask = and(theta > -pi/2, theta < pi/2);

subplot(1,3,2)
polarhistogram(theta,'BinEdges',linspace(0,2*pi,binnum))
% polarhistogram(theta_transform,'BinEdges',linspace(0,2*pi,binnum))
hold on
ax = gca;
polarplot(median(theta(fwd_mask),'omitnan')*[1,1],[0,ax.RLim(2)],'r--')
polarplot(median(theta(~fwd_mask),'omitnan')*[1,1],[0,ax.RLim(2)],'r--')
fprintf(strcat("Mean fwd theta (deg): ", num2str(mean(theta(fwd_mask)*180/pi,'omitnan')), "+/-", num2str(std(theta(fwd_mask)*180/pi,'omitnan')), "\n"))
set(gca,'ThetaZeroLocation','top')

% subplot(2,2,3)
% polarhistogram(theta(r>0),'BinWidth',bindeg*pi/180)
% set(gca,'ThetaZeroLocation','top')
% 
% subplot(2,2,4)
% polarhistogram(theta(r<0),'BinWidth',bindeg*pi/180)
% set(gca,'ThetaZeroLocation','top')

% Do a heat map on subplot 3

x = r.*cos(theta);
y = r.*sin(theta);


f2 = figure('Position',[300,300,600,1200]);
% subplot(1,3,3)
% hist3([y',x'], 'CDataMode','auto','Ctrs', {-20:5:20, -30:2:40});
hist3([y',x'], 'CDataMode','auto','Ctrs', {-20:5:20, -31.5:3:40.5});
% hist3([y',x'], 'CDataMode','auto','Ctrs', {-21:1.5:21, -33:1.5:42});
% hist3([y',x'], 'CDataMode','auto','Ctrs', {-21:3:21, -33:3:42});
% hist3([y',x'], 'CDataMode','auto','Ctrs', {-20:2.5:20, -33:3:42})
colormap('hot')
% caxis([0, max(max(z))]);

view(2)
% imagesc(flipud(hh3))
ax = gca;
set(ax,'XDir','reverse')
% ax.YTick = -40:8:50;
% ax.XTick = -20:8:20;
% set(ax, 'YTick',[0 yt], 'YTickLabel', [flip([0 yt])]*10)

end