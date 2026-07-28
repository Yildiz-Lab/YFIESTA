% https://www.mathworks.com/matlabcentral/answers/295595-arrange-multiple-2d-histograms-in-3d

%% Reading files, choose the AllStats file within each folder of the conditions you want to plot
% User change here what you want the tick labels to be
% list of figure paths to pull from


rootpath = '/Volumes/TOSHIBA EXT/Kif1a MAP7 Data/0nM Compiled Steps'; 
rootpath1 = '/Volumes/TOSHIBA EXT/Kif1a MAP7 Data/5nM Compiled Steps'; 
rootpath2 = '/Volumes/TOSHIBA EXT/Kif1a MAP7 Data/25nM Compiled Steps';
rootpath3 = '/Volumes/TOSHIBA EXT/Kif1a MAP7 Data/50nM Compiled Steps'; 
rootpath4 = '/Volumes/TOSHIBA EXT/Kif1a MAP7 Data/75nM Compiled Steps';
filename = 'AllStats.fig';
list_figures = {fullfile(rootpath,filename), fullfile(rootpath1,filename), fullfile(rootpath2,filename), fullfile(rootpath3,filename),fullfile(rootpath4,filename)};
x_tick_label = {'0','5','25','50','75'};
x_tick_label_arr = ["0","5","25","50","75"];


% %K490 MAP7 0819Leo
% rootpath = 'F:\Leo\Analyzed Stepping Data\K490 MAP7 Stepping 0812\0nM'; rootpath1 = 'F:\Leo\Analyzed Stepping Data\K490 MAP7 Stepping 0812\20nM'; rootpath2 = 'F:\Leo\Analyzed Stepping Data\K490 MAP7 Stepping 0812\50nM';
% rootpath3 = 'F:\Leo\Analyzed Stepping Data\K490 MAP7 Stepping 0812\100nM'; 
% filename = 'AllStats.fig';
% list_figures = {fullfile(rootpath,filename), fullfile(rootpath1,filename), fullfile(rootpath2,filename), fullfile(rootpath3,filename)};
% x_tick_label = {'0','6.6%','13%','17.6%'};
% x_tick_label_arr = ["0","6.6%","13%","17.6%"];


% % K560MAP9 0619Leo
% rootpath = 'E:\Leo\RAW data\20240613\1'; rootpath1 = 'E:\Leo\RAW data\20240613\2'; rootpath2 = 'E:\Leo\RAW data\20240613\3';
% rootpath3 = 'E:\Leo\RAW data\20240613\4'; rootpath4 = 'E:\Leo\RAW data\20240613\5'; rootpath5 = 'E:\Leo\RAW data\20240613\6';
% filename = 'AllStats.fig';
% list_figures = {fullfile(rootpath,filename), fullfile(rootpath1,filename), fullfile(rootpath2,filename), fullfile(rootpath3,filename),fullfile(rootpath4,filename),fullfile(rootpath5,filename)};
% x_tick_label = {'0','0.5%','1%','1.5%','2%','3%'};
% x_tick_label_arr = ["0","0.5%","1%","1.5%","2%","3%"];


% % K490MAP9 0819Leo
% rootpath = 'F:\Leo\Analyzed Stepping Data\K490 MAP9 Stepping 0814 0816\0nM'; rootpath1 = 'F:\Leo\Analyzed Stepping Data\K490 MAP9 Stepping 0814 0816\5nM'; rootpath2 = 'F:\Leo\Analyzed Stepping Data\K490 MAP9 Stepping 0814 0816\10nM';
% rootpath3 = 'F:\Leo\Analyzed Stepping Data\K490 MAP9 Stepping 0814 0816\20nM'; 
% filename = 'AllStats.fig';
% list_figures = {fullfile(rootpath,filename), fullfile(rootpath1,filename), fullfile(rootpath2,filename), fullfile(rootpath3,filename)};
% x_tick_label = {'0','0.5%','1%','2%'};
% x_tick_label_arr = ["0","0.5%","1%","2%"];


% % MINFLUX MAP7 20251209
% rootpath = 'E:\Leo\RAW data\20251119 MINFLUX\Summary\K490 MAP7\0nM';rootpath1 = 'E:\Leo\RAW data\20251119 MINFLUX\Summary\K490 MAP7\5nM'; 
% rootpath2 = 'E:\Leo\RAW data\20251119 MINFLUX\Summary\K490 MAP7\10nM'; 
% rootpath3 = 'E:\Leo\RAW data\20251119 MINFLUX\Summary\K490 MAP7\15nM';
% rootpath4 = 'E:\Leo\RAW data\20251119 MINFLUX\Summary\K490 MAP7\25nM'; rootpath5 = 'E:\Leo\RAW data\20251119 MINFLUX\Summary\K490 MAP7\50nM';
% filename = 'AllStats.fig';
% list_figures = {fullfile(rootpath,filename), fullfile(rootpath1,filename),fullfile(rootpath2,filename),fullfile(rootpath3,filename),fullfile(rootpath4,filename),fullfile(rootpath5,filename)};
% x_tick_label = {'0','9%','16.5%','22.9%','33.3%','50.5%'};
% x_tick_label_arr = ["0","9%","16.5%","22.9%","33.3%","50.5%"];
% x = 0.33*[0,0.2,0.5,0.71,0.86,0.93];

% % MINFLUX MAP7 20251209 without 50
% rootpath = 'F:\Leo\RAW data\20251119 MINFLUX\Summary\K490 MAP7\0nM';rootpath1 = 'F:\Leo\RAW data\20251119 MINFLUX\Summary\K490 MAP7\5nM'; 
% rootpath2 = 'F:\Leo\RAW data\20251119 MINFLUX\Summary\K490 MAP7\10nM'; 
% rootpath3 = 'F:\Leo\RAW data\20251119 MINFLUX\Summary\K490 MAP7\15nM';
% rootpath4 = 'F:\Leo\RAW data\20251119 MINFLUX\Summary\K490 MAP7\25nM'; 
% filename = 'AllStats.fig';
% list_figures = {fullfile(rootpath,filename), fullfile(rootpath1,filename),fullfile(rootpath2,filename),fullfile(rootpath3,filename),fullfile(rootpath4,filename)};
% x_tick_label = {'0','9%','16.5%','22.9%','33.3%'};
% x_tick_label_arr = ["0","9%","16.5%","22.9%","33.3%"];
% % x = 0.33*[0,0.2,0.5,0.71,0.86,0.93];

% % MINFLUX MAP9 20251209
% rootpath = 'F:\Leo\RAW data\20251119 MINFLUX\Summary\K490 MAP9\0nM';rootpath1 = 'F:\Leo\RAW data\20251119 MINFLUX\Summary\K490 MAP9\5nM'; 
% rootpath2 = 'F:\Leo\RAW data\20251119 MINFLUX\Summary\K490 MAP9\10nM'; 
% rootpath3 = 'F:\Leo\RAW data\20251119 MINFLUX\Summary\K490 MAP9\25nM'; rootpath4 = 'F:\Leo\RAW data\20251119 MINFLUX\Summary\K490 MAP9\40nM';
% filename = 'AllStats.fig';
% list_figures = {fullfile(rootpath,filename), fullfile(rootpath1,filename),fullfile(rootpath2,filename),fullfile(rootpath3,filename),fullfile(rootpath4,filename)};
% x_tick_label = {'0','2.2%','4.3%','10%','14.8%'};
% x_tick_label_arr = ["0","2.2%","4.3%","10%","14.8%"];
% % x = 0.33*[0,0.2,0.5,0.71];

% % MINFLUX K560 MAP7 20251209
% rootpath = 'E:\Leo\RAW data\20251119 MINFLUX\Summary\K560 MAP7\0nM';rootpath1 = 'E:\Leo\RAW data\20251119 MINFLUX\Summary\K560 MAP7\10nM'; 
% filename = 'AllStats.fig';
% list_figures = {fullfile(rootpath,filename), fullfile(rootpath1,filename)};
% x_tick_label = {'0','10nM'};
% x_tick_label_arr = ["0","10nM"];
% x = 0.33*[0,0.5];

% https://www.mathworks.com/matlabcentral/answers/295595-arrange-multiple-2d-histograms-in-3d
%% Reading files, choose the AllStats file within each folder of the conditions you want to plot
% User change here what you want the tick labels to be
% list of figure paths to pull from


%% Setup and Folder Selection
% Ask user where to save the output figures
save_dir = uigetdir(pwd, 'Select a folder to save output figures');
if save_dir == 0
    error('No folder selected. Script aborted.');
end

on_axis = cell(1,length(list_figures));
off_axis = cell(1,length(list_figures));
for i = 1:length(list_figures)
    f = open(list_figures{i});
    har = get(f,'children');
    %har{9} should be on-axis steps,s har{7} off
    on_axis{i} = har(9).Children.Data;
    off_axis{i} = har(7).Children.Data;
    close(f)
end

%% Adapting and packaging
%parameters for plotting, binedges, etc.
on_binedges = -40:3:40;
off_binedges = -40:3:40;
sb = 7; %space between consecutive histograms filled with white space
tl = sb*length(list_figures)-sb+1;
z_on = zeros(length(on_binedges)-1, tl);
z_off = zeros(length(off_binedges)-1, tl);
p_on = zeros(1,length(list_figures));
p_off = zeros(1,length(list_figures));
ss = zeros(2,length(list_figures));
large_steps_on = zeros(1,length(list_figures));
large_steps_off = zeros(1,length(list_figures));
ms_on = zeros(1,length(list_figures));
c1s_on = zeros(1,length(list_figures));
c2s_on = zeros(1,length(list_figures));
ms_off = zeros(1,length(list_figures));
c1s_off = zeros(1,length(list_figures));
c2s_off = zeros(1,length(list_figures));

for i = 1:length(list_figures)
    h = histogram(on_axis{i},on_binedges);
    z_on(:,sb*i-sb+1) = h.Values';
    total_parallel_steps = length(h.Data);
    total_fwd_parallel_steps = sum(h.Data > 0);
    large_steps_on(i) = sum(abs(h.Data) > 32)/length(h.Data);
    h = histogram(off_axis{i},off_binedges);
    z_off(:,sb*i-sb+1) = h.Values';
    total_transverse_steps = length(h.Data);
    large_steps_off(i) = sum(abs(h.Data) > 16)/total_transverse_steps;
    
    % Normalization
    z_on(:,sb*i-sb+1) = z_on(:,sb*i-sb+1)/total_fwd_parallel_steps;
    z_off(:,sb*i-sb+1) = z_off(:,sb*i-sb+1)/total_fwd_parallel_steps;
    
    % percentage back or sideways
    Nback = total_parallel_steps - total_fwd_parallel_steps;
    p_on(i) = Nback / total_fwd_parallel_steps;
    p_off(i) = length(off_axis{i}) / total_fwd_parallel_steps;

    [ms_on(i),c1s_on(i),c2s_on(i)] = beta_confidence(Nback,total_parallel_steps+total_transverse_steps-Nback);
    [ms_off(i),c1s_off(i),c2s_off(i)] = beta_confidence(total_transverse_steps,total_parallel_steps);
    mdl_gauss = fittype('normcdf(x,mu,sigma)','indep','x');
    X = sort(abs(h.Data));
    Y = linspace(0,1,length(X));
    fittedmdl = fit(X, Y', mdl_gauss, ...
    'StartPoint', [8, 8], ...
    'Lower', [-Inf, 0.01], ...  % sigma >= 0.01
    'Upper', [Inf, Inf]);
    ss(1,i) = fittedmdl.mu; ss(2,i) = fittedmdl.sigma;
end

% Initialize error arrays
sem_on = zeros(1,length(list_figures));
sem_off = zeros(1,length(list_figures));

%error bars
for i = 1:length(list_figures)
    h = histogram(on_axis{i},on_binedges);
    z_on(:,sb*i-sb+1) = h.Values';
    total_parallel_steps = length(h.Data);
    total_fwd_parallel_steps = sum(h.Data > 0);
    
    h_off = histogram(off_axis{i},off_binedges);
    z_off(:,sb*i-sb+1) = h_off.Values';
    total_transverse_steps = length(h_off.Data);
    
    % --- NORMALIZATION ---
    z_on(:,sb*i-sb+1) = z_on(:,sb*i-sb+1)/total_fwd_parallel_steps;
    z_off(:,sb*i-sb+1) = z_off(:,sb*i-sb+1)/total_fwd_parallel_steps;
    
    % --- BAR VALUES ---
    Nback = total_parallel_steps - total_fwd_parallel_steps;
    p_on(i) = Nback / total_fwd_parallel_steps;
    p_off(i) = total_transverse_steps / total_fwd_parallel_steps;
    
    % --- ERROR BARS ---
    if Nback > 0
        sem_on(i) = p_on(i) * sqrt(1/Nback + 1/total_fwd_parallel_steps);
    else
        sem_on(i) = 0; 
    end
    
    if total_transverse_steps > 0
        sem_off(i) = p_off(i) * sqrt(1/total_transverse_steps + 1/total_fwd_parallel_steps);
    else
        sem_off(i) = 0;
    end
end

%% Plotting - 3D On-Axis
fig1 = figure();
z = z_on;
h = bar3(on_binedges(1:end-1),z);
view([-45 25]);

% BOLD AND ENLARGE
title('On-axis Step Distribution...','FontSize',16,'FontWeight','bold')
xlabel('Density(%)','FontSize',14,'FontWeight','bold')              
ylabel('On-axis Step (nm)','FontSize',14,'FontWeight','bold')       
zlabel('Normalized Frequency','FontSize',14,'FontWeight','bold')    
set(gca,'FontSize',12,'LineWidth',1.5,'FontWeight','bold')          

cm_choice = colormap(winter);
interpolation = 30; 
cm = ones(tl,3);
ecm = ones(tl,3);
for i = 1:tl
    if mod(i-1,sb) == 0 
        ii = ceil(i/sb);
        cm(i,:) = cm_choice(interpolation*(ii+1),:);
        ecm(i,:) = zeros(1,3);
    end
end
set(gca,'XTick',1:sb:tl)
set(gca,'XTickLabel',x_tick_label)
set(gca,'YLim',[min(on_binedges)-sb,max(on_binedges)])
set(gca,'YTick',16*sign(min(on_binedges))*floor(abs(min(on_binedges)/16)):16:16*floor(max(on_binedges)/16))
set(gca,'ZLim',[0,2*max(max([h.ZData]))]) 

cnt = 0;
for jj = 1:length(h)
    xd = get(h(jj),'xdata');
    yd = get(h(jj),'ydata');
    zd = get(h(jj),'zdata');
    delete(h(jj))    
    idx = [0;find(all(isnan(xd),2))];
    if jj == 1
        S = zeros(length(h)*(length(idx)-1),1);
        dv = floor(size(cm,1)/length(S));
    end
    for ii = 1:length(idx)-1
        cnt = cnt + 1;
        S(cnt) = surface(xd(idx(ii)+1:idx(ii+1)-1,:),...
                         yd(idx(ii)+1:idx(ii+1)-1,:),...
                         zd(idx(ii)+1:idx(ii+1)-1,:),...
                         'facecolor',cm(ceil(cnt/size(z,1)),:),...
                         'EdgeColor',ecm(ceil(cnt/size(z,1)),:),...
                         'FaceAlpha',0.9);
    end
end
rotate3d
grid off

% SAVE FIGURE 1
saveas(fig1, fullfile(save_dir, 'OnAxis_3D_Histogram.fig')); % Keep .fig as is
exportgraphics(fig1, fullfile(save_dir, 'OnAxis_3D_Histogram.pdf'), 'Resolution', 300); % Change PDF line
fprintf('Saved OnAxis_3D_Histogram\n');


%% Plotting - 3D Off-axis
fig2 = figure();
z = z_off;
h = bar3(off_binedges(1:end-1),z);
view([-45 25]);

% BOLD AND ENLARGE
title('Off-axis Step Distribution (3D Histogram)','FontSize',16,'FontWeight','bold')
xlabel('Density(%)','FontSize',14,'FontWeight','bold')
ylabel('Off-axis Step (nm)','FontSize',14,'FontWeight','bold')
zlabel('Normalized Frequency','FontSize',14,'FontWeight','bold')
set(gca,'FontSize',12,'LineWidth',1.5,'FontWeight','bold') 

set(gca,'XTick',1:sb:tl)
set(gca,'XTickLabel',x_tick_label)
set(gca,'YLim',[min(off_binedges)-sb,max(off_binedges)])
set(gca,'YTick',16*sign(min(off_binedges))*floor(abs(min(off_binedges)/16)):16:16*floor(max(off_binedges)/16))
set(gca,'ZLim',[0,2.5*max(max([h.ZData]))]) 

cnt = 0;
for jj = 1:length(h)
    xd = get(h(jj),'xdata');
    yd = get(h(jj),'ydata');
    zd = get(h(jj),'zdata');
    delete(h(jj))    
    idx = [0;find(all(isnan(xd),2))];
    if jj == 1
        S = zeros(length(h)*(length(idx)-1),1);
        dv = floor(size(cm,1)/length(S));
    end
    for ii = 1:length(idx)-1
        cnt = cnt + 1;
        S(cnt) = surface(xd(idx(ii)+1:idx(ii+1)-1,:),...
                         yd(idx(ii)+1:idx(ii+1)-1,:),...
                         zd(idx(ii)+1:idx(ii+1)-1,:),...
                         'facecolor',cm(ceil(cnt/size(z,1)),:),...
                         'EdgeColor',ecm(ceil(cnt/size(z,1)),:),...
                         'FaceAlpha',0.9);
    end
end
rotate3d
grid off

% SAVE FIGURE 2
saveas(fig2, fullfile(save_dir, 'OffAxis_3D_Histogram.fig'));
exportgraphics(fig2, fullfile(save_dir, 'OffAxis_3D_Histogram.pdf'), 'Resolution', 300);
fprintf('Saved OffAxis_3D_Histogram\n');


%% Plots for side & backward step percentages bar plots
X = 1:length(list_figures);

% --- Fraction Backward ---
fig3 = figure();
hold on
c_on = bar(p_on);

% BOLD AND ENLARGE
title('Fraction of Backward On-axis Steps','FontSize',18,'FontWeight','bold')
xlabel('Density(%)','FontSize',16,'FontWeight','bold')
ylabel('Fraction','FontSize',16,'FontWeight','bold')
set(gca,'FontSize',14,'LineWidth',1.5,'FontWeight','bold')

c_on.FaceColor = 'flat';
for i = 1:length(p_on)
    c_on.CData(i,:) = cm_choice(interpolation*(i+1),:);
end
errorbar(X, p_on, sem_on, sem_on, 'k', 'LineStyle', 'none', 'LineWidth', 1.8, 'CapSize', 10);

% Dynamic Ylim
% max_val = max(p_on + sem_on);
% 
% % Set limit to 1.2x the max value (20% padding). 
% % If max is 0 (empty data), default to 0.1 so graph isn't flat.
% if max_val == 0
%     upper_lim = 0.1; 
% else
%     upper_lim = max_val * 1.2;
% end
% set(gca,'YLim',[0, upper_lim]) 
% % ------------------------------

set(gca,'YLim',[0,0.5])
set(gca, 'XTick', 1:length(p_on), 'XTickLabel', x_tick_label); 

% SAVE FIGURE 3
saveas(fig3, fullfile(save_dir, 'Fraction_Backward_Steps.fig'));
exportgraphics(fig3, fullfile(save_dir, 'Fraction_Backward_Steps.pdf'), 'Resolution', 300);
exportgraphics(fig3, fullfile(save_dir, 'Fraction_Backward_Steps.png'), 'Resolution', 300);
fprintf('Saved Fraction_Backward_Steps\n');


% --- Fraction Off-axis ---
fig4 = figure();
hold on
c_off = bar(p_off);

% BOLD AND ENLARGE
title('Fraction of Off-axis Steps','FontSize',18,'FontWeight','bold')
xlabel('Density(%)','FontSize',16,'FontWeight','bold')
ylabel('Fraction','FontSize',16,'FontWeight','bold')
set(gca,'FontSize',14,'LineWidth',1.5,'FontWeight','bold')

c_off.FaceColor = 'flat';
for i = 1:length(p_on)
    c_off.CData(i,:) = cm_choice(interpolation*(i+1),:);
end
errorbar(X, p_off, sem_off, sem_off, 'k', 'LineStyle', 'none', 'LineWidth', 1.8, 'CapSize', 10);

% For dynamic Ylim
max_val = max(p_off + sem_off);

% Set limit to 1.2x the max value (20% padding)
if max_val == 0
    upper_lim = 0.1;
else
    upper_lim = max_val * 1.2;
end
set(gca,'YLim',[0, upper_lim])
% ------------------------------
% set(gca,'YLim',[0,0.6])
set(gca, 'XTick', 1:length(p_off), 'XTickLabel', x_tick_label_arr);

% SAVE FIGURE 4
saveas(fig4, fullfile(save_dir, 'Fraction_OffAxis_Steps.fig'));
exportgraphics(fig4, fullfile(save_dir, 'Fraction_OffAxis_Steps.pdf'), 'Resolution', 300);
exportgraphics(fig4, fullfile(save_dir, 'Fraction_OffAxis_Steps.png'), 'Resolution', 300);
fprintf('Saved Fraction_OffAxis_Steps\n');


%% Plots for side step size
fig5 = figure();
hold on
c_on = bar(1:length(p_on),ss(1,:));

% BOLD AND ENLARGE
title('Fitted Step Size (CDF Gaussian Fit)','FontSize',18,'FontWeight','bold')
xlabel('Density(%)','FontSize',16,'FontWeight','bold')
ylabel('Step Size (nm)','FontSize',16,'FontWeight','bold')
set(gca,'FontSize',14,'LineWidth',1.5,'FontWeight','bold')

c_on.FaceColor = 'flat';
for i = 1:length(p_on)
    c_on.CData(i,:) = cm_choice(interpolation*(i+1),:);
end
er = errorbar(1:length(p_on),ss(1,:),ss(2,:),'Color',[0 0 0],'LineStyle','none');
set(gca, 'XTick', 1:length(p_on), 'XTickLabel', x_tick_label_arr);

% SAVE FIGURE 5
saveas(fig5, fullfile(save_dir, 'Step_Size_Fitted.fig'));
exportgraphics(fig5, fullfile(save_dir, 'Step_Size_Fitted.pdf'), 'Resolution', 300);
fprintf('Saved Step_Size_Fitted\n');

fprintf('All figures saved successfully to: %s\n', save_dir);