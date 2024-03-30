function [nearestPoints] = findNearestPointOnLine_chatGPT(x_line, y_line)

% Plot the line
    f = figure('KeyPressFcn', @keyPressCallback);
    ax = axes('Parent', f);

    % Initialize for storage and selecting points
    f.UserData.nearestPoints = [];
    f.UserData.selectPoints = false;
    
    plot(x_line, y_line, 'b.-', 'tag', 'dominant axis');
    hold on;  
    
    uiwait(f)

% end

% Callback function for key press
function keyPressCallback(f, event)
    if strcmp(event.Key, 'return')
        f.UserData.selectPoints = ~f.UserData.selectPoints; % Toggle selection flag
        if f.UserData.selectPoints
            disp('Selecting points: ON');
            hold on;
            % Continue until Enter key is pressed
            while true
                % Prompt user to click
                title('Click on the plot to find the nearest point on the line. Press Enter to stop.');
                [x_click, y_click, key] = ginput(1);

                % Check if Enter key is pressed
                if strcmp(key, 'return')
                    disp('return \n')
                    break; %return
                end

                try
                    % Plot the clicked point
%                     plot(x_click, y_click, 'ro');
                    
                    ax = f.Children;
                    dominant_line = findobj(ax,'tag','dominant axis');
                    x_line = dominant_line.XData;
                    y_line = dominant_line.YData;

                    % Find the nearest point on the line
                    distances = sqrt((x_line - x_click).^2);
                    [~, index] = min(distances);
                    nearestPoint = [x_line(index), y_line(index)];

                    % Plot the nearest point
                    plot(nearestPoint(1), nearestPoint(2), 'go', 'MarkerSize', 10);

                    % Store the nearest point for calculation purposes
                    f.UserData.nearestPoints = [f.UserData.nearestPoints; nearestPoint];
                    nearestPoint
                    
                    
                catch
                    break;
                end

            end
            
        else
            title('Not selecting points.');
            disp('Selecting points: OFF');
        end
    end
end

end


