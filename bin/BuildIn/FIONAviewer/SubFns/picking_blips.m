function picking_blips

% JS quick function 2024/03/20 to load previously analyzed fiona files and go back
% through to pick out blips. Save to data.blip the timepoint where it
% occurs for later analysis of average change in time and distance.

% This is freaking ugly, don't tell mom

[savefile, savefolder] = uigetfile('*.mat');
savename = fullfile(savefolder,savefile);
fdata = load(savename);
[pathname, filename, ext] = fileparts(savename);
savename = fullfile(pathname,strcat(filename,'_wblip',ext));

data = fdata.data;

x_line = data.time;
y_line = data.trace(:,1);

% Plot the line
f = figure('KeyPressFcn', @keyPressCallback);
ax = axes('Parent', f);


hold on;
plot(x_line, y_line, 'b.-', 'tag', 'dominant axis');
plot(x_line, data.trace(:,3), 'g-', 'tag', 'dominant fit');

% calculate std and plot 2std
sigma = std(data.trace(:,1) - data.trace(:,3),'omitnan')
plot(x_line, data.trace(:,3)+2*sigma, 'k--', 'tag', '2sigma')
plot(x_line, data.trace(:,3)-2*sigma, 'k--', 'tag', '2sigma')

% Initialize for storage and selecting points
f.UserData.selectPoints = false;
f.UserData.sigma = sigma;
if isfield(data,'blips')
    f.UserData.nearestPoints = data.blips;
    % blips already started, appending
    plot(f.UserData.nearestPoints(:,2), f.UserData.nearestPoints(:,3), 'go', 'MarkerSize', 10);
else
    f.UserData.nearestPoints = [];
    data.blips = [];
end

save(savename, 'data')

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
                    nearestPoint = [index, x_line(index), y_line(index)];

                    % Plot the nearest point
                    plot(nearestPoint(2), nearestPoint(3), 'go', 'MarkerSize', 10);

                    % Store the nearest point for calculation purposes
                    f.UserData.nearestPoints = [f.UserData.nearestPoints; nearestPoint];
                    nearestPoint
                    data.blips = f.UserData.nearestPoints;
                    save(savename, 'data')
                    
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


