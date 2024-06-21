function plot_step_mod_mean(fwd_mod_mean, bwd_mod_mean, time_bins)
%fwd and bwd_mod_mean should be cells (pertaining to traces) of (Nx2) cells
%(pertaining to specific steps) where the first column is modded time
%aligned to step and the second column is modded with average step position

if nargin < 3
    time_bins = 0:0.0005:0.05;
end

% first, we should extract the cells into one large cell, so:
fwd_extracted = cell(0,2);
bwd_extracted = cell(0,2);
for j = 1:length(fwd_mod_mean)
    trace_steps = fwd_mod_mean{j};
    m = length(fwd_extracted);
    for i = 1:length(trace_steps)
        fwd_extracted{m+i,1} = trace_steps{i,1};
        fwd_extracted{m+i,2} = trace_steps{i,2};
    end
    trace_steps = bwd_mod_mean{j};
    m = length(bwd_extracted);
    for i = 1:length(trace_steps)
        bwd_extracted{m+i,1} = trace_steps{i,1};
        bwd_extracted{m+i,2} = trace_steps{i,2};
    end
end

% Now we have all time in column 1 and position in column 2
% extract and plot them  all over each other.
fwd_array = zeros(0,2);
figure()
hold on
for j = 1:length(fwd_extracted)
    plot(fwd_extracted{j,1},fwd_extracted{j,2});
    fwd_array = [fwd_array; fwd_extracted{j,1}, fwd_extracted{j,2}];
end
% for j = 1:length(bwd_extracted)
%     plot(bwd_extracted{j,1},bwd_extracted{j,2});
%     bwd_array = [bwd_array; bwd_extracted{j,1}, bwd_extracted{j,2}];
% end

meanpos = zeros(1,length(time_bins)-1);
stdpos = zeros(1,length(time_bins)-1);
for k = 1:length(time_bins)-1
    timemask = and(fwd_array(:,1) > -time_bins(k+1), fwd_array(:,1) < -time_bins(k));
    meanpos(k) = mean(fwd_array(timemask,2),'omitnan');
    stdpos(k) = std(fwd_array(timemask,2),'omitnan');
end

t = -time_bins(1:end-1)-(time_bins(2)-time_bins(1))/2;
plot(t,meanpos,'k-','LineWidth',5)
plot(t,meanpos+stdpos,'Color',[0.3,0.3,0.3],'LineWidth',2)
plot(t,meanpos-stdpos,'Color',[0.3,0.3,0.3],'LineWidth',2)
fill([t, fliplr(t)], [meanpos+stdpos, fliplr(meanpos-stdpos)],[0.5,0.5,0.5],'EdgeColor','None','FaceAlpha',0.3)

text(-time_bins(round(end/2)),30,strcat("N = ",num2str(length(fwd_extracted))))
ax = gca;
ax.XLim = sort([-min(time_bins)+0.005,-max(time_bins)-0.005]);


bwd_array = zeros(0,2);
figure()
hold on
for j = 1:length(bwd_extracted)
    plot(bwd_extracted{j,1},bwd_extracted{j,2});
    bwd_array = [bwd_array; bwd_extracted{j,1}, bwd_extracted{j,2}];
end

meanpos = zeros(1,length(time_bins)-1);
stdpos = zeros(1,length(time_bins)-1);
for k = 1:length(time_bins)-1
    timemask = and(bwd_array(:,1) > time_bins(k), bwd_array(:,1) < time_bins(k+1));
    meanpos(k) = mean(bwd_array(timemask,2),'omitnan');
    stdpos(k) = std(bwd_array(timemask,2),'omitnan');
end

plot(time_bins(1:end-1)+(time_bins(2)-time_bins(1))/2,meanpos,'k-','LineWidth',5)
plot(time_bins(1:end-1)+(time_bins(2)-time_bins(1))/2,meanpos+stdpos,'k--','LineWidth',2)
plot(time_bins(1:end-1)+(time_bins(2)-time_bins(1))/2,meanpos-stdpos,'k--','LineWidth',2)

text(time_bins(round(end/2)),30,strcat("N = ",num2str(length(bwd_extracted))))
ax = gca;
ax.XLim = sort([min(time_bins)-0.005,max(time_bins)+0.005]);

end