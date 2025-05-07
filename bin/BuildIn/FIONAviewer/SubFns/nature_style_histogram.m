function nature_style_histogram(data)

if nargin < 1
% Sample data
data = randn(1, 1000); % Example: normally distributed data
end

% Create histogram
figure;
h = histogram(data, 'FaceColor', 'k', 'EdgeColor', 'k');

% Set axis properties
set(gca, ...
    'FontName', 'Arial', ...           % Arial font
    'FontSize', 10, ...                % Reasonable size for Nature
    'TickDir', 'out', ...              % Ticks outside
    'LineWidth', 1, ...                % Thicker axis lines
    'Box', 'off', ...                  % Remove top/right box lines
    'XColor', 'k', ...                 % X axis black
    'YColor', 'k');                    % Y axis black

% Labels (customize if needed)
xlabel('X axis label', 'FontName', 'Arial', 'FontSize', 10, 'Color', 'k');
ylabel('Count', 'FontName', 'Arial', 'FontSize', 10, 'Color', 'k');

% Optional: Adjust figure background
set(gcf, 'Color', 'w'); % White background


end