% Get current axes
ax = gca;

% Get current x limits
xl = xlim(ax);
yl = ylim(ax);

% Find all line objects in the axes
hLines = findobj(ax, 'Type', 'line');

% Create new figure
figure(4);
hold on;
% grid on;

for k = 1:length(hLines)
    x = get(hLines(k),'XData');
    y = get(hLines(k),'YData');

    % Clip to visible X-range
    mask = x >= xl(1) & x <= xl(2);
    x_clip = x(mask);
    y_clip = y(mask);

    if isempty(x_clip)
        continue;
    end
    
    % Shift so origin is (0,0)
    x_shift = x_clip - xl(1);      % shift x to 0
    y_shift = y_clip - yl(1);      % shift y to 0

    % Plot with same style
    plot(x_shift, y_shift, ...
        'Color', get(hLines(k),'Color'), ...
        'LineStyle', get(hLines(k),'LineStyle'), ...
        'LineWidth', get(hLines(k),'LineWidth'), ...
        'Marker', get(hLines(k),'Marker'));
end

box off
set(gca, ...
    'FontName','Arial', ...
    'FontSize', 10, ...
    'LineWidth', 0.75, ...
    'TickDir','out', ...
    'TickLength',[0.015 0.015])

xlabel('Shifted Time (s)');
ylabel('Shifted Position');
title('Replotted Lines (Shifted to Origin)');
hold off;
