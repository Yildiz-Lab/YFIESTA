% https://www.mathworks.com/matlabcentral/answers/295595-arrange-multiple-2d-histograms-in-3d

%% Reading
rootpath = '/Volumes/Ultra Touch/231217';
filename = 'AllStats.fig';

list_figures = {fullfile(rootpath,'0nM RK',filename),fullfile(rootpath,'0.8nM RK',filename), fullfile(rootpath,'2nM RK',filename), fullfile(rootpath,'5nM RK',filename)};
on_axis = cell(1,length(list_figures));
off_axis = cell(1,length(list_figures));

for i = 1:length(list_figures)
    f = open(list_figures{i});
    har = get(f,'children');
    %har{9} should be on-axis steps, har{7} off
    on_axis{i} = har(9).Children.Data;
    off_axis{i} = har(7).Children.Data;
    close(f)
end

%% Adapting and packaging

%parameters for plotting, binedges, etc.
on_binedges = -24:2:66;
off_binedges = -42:2:42;
sb = 5; %space between consecutive histograms

% the 2x allows for a space between, visually more appealing
tl = sb*length(list_figures)-sb+1;
z_on = zeros(length(on_binedges)-1, tl);
z_off = zeros(length(off_binedges)-1, tl);

for i = 1:length(list_figures)
    h = histogram(on_axis{i},on_binedges);
    z_on(:,sb*i-sb+1) = h.Values';
    h = histogram(off_axis{i},off_binedges);
    z_off(:,sb*i-sb+1) = h.Values';
    % option to normalize
    total_parallel_steps = sum(z_on(:,sb*i-sb+1));
    z_on(:,sb*i-sb+1) = z_on(:,sb*i-sb+1)/total_parallel_steps;
    z_off(:,sb*i-sb+1) = z_off(:,sb*i-sb+1)/total_parallel_steps;
end

% for i = 1:length(list_figures)
%     hist(on_axis{i});
%     h = get(gca,'Children'); 
%     x = h.Vertices(:,1);
%     y = h.Vertices(:,2);
%     z = 3*i*ones(size(x));
%     close all
%     patch(x,y,z, 'r'); hold on;
% end
% view(3)



%% Plotting

% z = [1 0 4 0 7; 2 0 5 0 8; 3 0 6 0 9; 4 0 7 0 10];
z = z_on;
h = bar3(on_binedges(1:end-1),z);
x_tick_label = {'0','0.8 [0.045]','2 [0.11]','5 [0.23]'};

x_ticks_mod = cell(1,tl);

%set(gca,'YTickLabel',on_binedges-2)
% cm = get(gcf,'colormap');  % Use the current colormap.
cm_choice = colormap(winter);
interpolation = 60; %round(length(cm)/length(list_figures));
% generate colormap
cm = ones(tl,3);
ecm = ones(tl,3);
for i = 1:tl
    if mod(i-1,sb) == 0 %else fill in blank space
        ii = ceil(i/sb);
        cm(i,:) = cm_choice(interpolation*ii,:);
        ecm(i,:) = zeros(1,3);
        x_ticks_mod{i} = x_tick_label{ii};
    end
end

set(gca,'XTickLabel',x_ticks_mod)
% cm = [0.5 0.5 0.5; 1 1 1; 0.5 0.5 0.5; 1 1 1; 0.5 0.5 0.5];
% ecm = zeros(tl,3); ecm(2:2:tl,:) = ones(length(2:2:tl),3);

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
    % for ii = 1:length(idx)-1
    %     cnt = cnt + 1;
    %     S(cnt) = surface(xd(idx(ii)+1:idx(ii+1)-1,:),...
    %                      yd(idx(ii)+1:idx(ii+1)-1,:),...
    %                      zd(idx(ii)+1:idx(ii+1)-1,:),...
    %                      'facecolor',cm((cnt-1)*dv+1,:));
    % end
    % for ii = 1:length(idx)-1
    %     cnt = cnt + 1;
    %     S(cnt) = surface(xd(idx(ii)+1:idx(ii+1)-1,:),...
    %                      yd(idx(ii)+1:idx(ii+1)-1,:),...
    %                      zd(idx(ii)+1:idx(ii+1)-1,:),...
    %                      'facecolor',cm(ceil(cnt/size(z,1)),:),...
    %                      'EdgeColor',ecm(ceil(cnt/size(z,1)),:));
    % end
    for ii = 1:length(idx)-1
        cnt = cnt + 1;
        S(cnt) = surface(xd(idx(ii)+1:idx(ii+1)-1,:),...
                         yd(idx(ii)+1:idx(ii+1)-1,:),...
                         zd(idx(ii)+1:idx(ii+1)-1,:),...
                         'facecolor',cm(ceil(cnt/size(z,1)),:),...
                         'EdgeColor',ecm(ceil(cnt/size(z,1)),:),...
                         'FaceAlpha',0.7);
    end
end
rotate3d

% now off_axis
figure()
z = z_off;
h = bar3(off_binedges(1:end-1),z);
set(gca,'XTickLabel',x_ticks_mod)

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
    % for ii = 1:length(idx)-1
    %     cnt = cnt + 1;
    %     S(cnt) = surface(xd(idx(ii)+1:idx(ii+1)-1,:),...
    %                      yd(idx(ii)+1:idx(ii+1)-1,:),...
    %                      zd(idx(ii)+1:idx(ii+1)-1,:),...
    %                      'facecolor',cm((cnt-1)*dv+1,:));
    % end
    % for ii = 1:length(idx)-1
    %     cnt = cnt + 1;
    %     S(cnt) = surface(xd(idx(ii)+1:idx(ii+1)-1,:),...
    %                      yd(idx(ii)+1:idx(ii+1)-1,:),...
    %                      zd(idx(ii)+1:idx(ii+1)-1,:),...
    %                      'facecolor',cm(ceil(cnt/size(z,1)),:),...
    %                      'EdgeColor',ecm(ceil(cnt/size(z,1)),:));
    % end
    for ii = 1:length(idx)-1
        cnt = cnt + 1;
        S(cnt) = surface(xd(idx(ii)+1:idx(ii+1)-1,:),...
                         yd(idx(ii)+1:idx(ii+1)-1,:),...
                         zd(idx(ii)+1:idx(ii+1)-1,:),...
                         'facecolor',cm(ceil(cnt/size(z,1)),:),...
                         'EdgeColor',ecm(ceil(cnt/size(z,1)),:),...
                         'FaceAlpha',0.7);
    end
end
rotate3d


%% Plots for side & backward step percentages bar plots



