% plotting random stuff for two color, eventually wrap into FIESTA

% load('/Volumes/TOSHIBA EXT/MINFLUX JS/241109/241109_03/241109-151001_minflux_Joseph_interleaved_Cy3controller_bg75kCy3Cy5_fiesta.mat')
% Molecule6293 = Molecule(78).Results;
% Molecule6295 = Molecule(79).Results;
% 
% 
% % Define the range of time
% time = linspace(0, 10, 100);  % Time points
% 
% % Define the x, y, and z (time) coordinates
% x6293 = Molecule6293(4:end,3); x6295 = Molecule6295(4:end,3);              % Example: x varies as sin(time)
% y6293 = Molecule6293(4:end,4); y6295 = Molecule6295(4:end,4);              % Example: y varies as cos(time)
% z6293 = Molecule6293(4:end,2); z6295 = Molecule6295(4:end,2);                   % z-coordinate is time
% 
% % Plot the 3D line
% figure;
% plot3(x6293, y6293, z6293, 'Color', 'g', 'LineWidth', 2);  % The 'LineWidth' argument makes the line thicker
% hold on
% plot3(x6295, y6295, z6295, 'Color', [1 0 1], 'LineWidth', 2);
% 
% % Add labels and title
% xlabel('X');
% ylabel('Y');
% zlabel('Time');
% title('3D Line Plot of Molecule 6293 & 6295 (X, Y, Time)');
% 
% % Adjust the view angle for better visualization
% view(45, 30);
% grid on;  % Add grid for better visualization
% 
% 
% figure;
% plot3(x6293, y6293, z6293, 'Color', 'g', 'LineWidth', 2, 'LineStyle', '--');  % The 'LineWidth' argument makes the line thicker
% scatter3(x6293, y6293, z6293, 25, 'g', 'filled')
% hold on
% plot3(x6295, y6295, z6295, 'Color', [1 0 1], 'LineWidth', 2, 'LineStyle', '--');  % The 'LineWidth' argument makes the line thicker
% scatter3(x6295, y6295, z6295, 25, [1 0 1], 'filled')
% 
% % Add labels and title
% xlabel('X');
% ylabel('Y');
% zlabel('Time');
% title('3D Line Plot of Molecule 6293 & 6295 (X, Y, Time)');
% 
% % Adjust the view angle for better visualization
% view(45, 30);
% grid on;  % Add grid for better visualization


%% The next one 50k 100k K490

% load('/Volumes/TOSHIBA EXT/MINFLUX JS/241109/241109_03/241109-151830_minflux_Joseph_interleaved_Cy3controller_bg50kCy3_100kCy5_fiesta.mat')
load('/Volumes/TOSHIBA EXT/MINFLUX JS/241109/241109_03/241109-151830_bg50kCy3_100kCy5_fiesta_path.mat')
load('/Users/slivka/Documents/MATLAB/241109-151830_minflux_Joseph_interleaved_Cy3controller_bg50kCy3_100kCy5_fiesta_corrected.mat')
load('/Users/slivka/Documents/MATLAB/241109-151830_minflux_fiesta_nocorrect.mat')
load('/Users/slivka/Documents/MATLAB/241109-151830_minflux_fiesta_pluscorrect.mat')

% 28612, 28613 -> 567, 568 [-7, -2]
% 26472, 26474 -> 523, 524 
% 24741, 24742 -> 487, 488
% 24445, 24446 -> 473, 474
% 22308, 22310 -> 439, 440
% 21993, 21994 -> 421, 422
% 14679, 14682 -> 289, 290
% 8848, 8842 -> 175, 176

% %% Or some dynein
% load('/Volumes/TOSHIBA EXT/MINFLUX JS/241112/241112_03/241112-142335_minflux_655_DPT_555_DPJ_LP2p5_fiesta_wpath.mat')
% 
% % 156, 155 -> 4,3
% % 2553, 2552 -> 20,19
% % 10555, 10551 -> 90,88
% % 33336, 33335 -> 274, 273
% % 68525, 68524 -> 493, 492
% 
% load('/Volumes/TOSHIBA EXT/MINFLUX JS/241112/241112_03/241112-141511_minflux_555_DPJ_655_DPT_LP2p5_fiesta_wpath2.mat')
% % 1110, 1112 -> 6, 7
% % 25542, 25543 -> 193, 194
% 
% load('/Volumes/TOSHIBA EXT/MINFLUX JS/241112/241112_05/241112-173933_minflux_555_JSfast_LP3_655_JSfast_LP8_fiesta_wpath.mat')
% % 4815, 4816 -> 21, 22
% % 16937, 16942 -> 49, 50
% % 28048, 28049 -> 71, 72
% % 77830, 77847 -> 173, 174
% % 83818, 83819 -> 186, 187
% % 104493, 104494 -> 228, 229
% 
% load('/Volumes/TOSHIBA EXT/MINFLUX JS/241112/241112_05/241112-174726_minflux_555_JSfast_LP4_655_JSfast_LP8_fiesta_wpath.mat')
% % 54368, 54371 -> 107, 108
% % 70934, 70935 -> 134, 135
% 
% load('/Volumes/TOSHIBA EXT/MINFLUX JS/241112/241112_05/241112-175641_minflux_555_JSfast_LP4_655_JSfast_LP8_fiesta_wpath.mat')
% % 30844, 30845 -> 63, 64
% 
% load('/Volumes/TOSHIBA EXT/MINFLUX JS/241112/241112_05/241112-180557_minflux_555_JSfast_LP4_655_JSfast_LP9_fiesta_wpath.mat')
% % 10087, 10088 -> 16, 17
% % 48357, 48358 -> 93, 94
% 
% load('/Volumes/TOSHIBA EXT/MINFLUX JS/241113/241113_01/241113-100126_minflux_fiesta_wpath.mat')
% % 47094, 47095 -> 90, 91
% 
% load('/Volumes/TOSHIBA EXT/MINFLUX JS/241113/241113_03/241113-112331_minflux_dcop_fiesta_wpath.mat')
% % 12540, 12542 -> 39, 40
% % 12933, 12934 -> 41, 42
% % 53135, 53137 -> 144, 145
% 
% load('/Volumes/TOSHIBA EXT/MINFLUX JS/241113/241113_03/241113-112810_minflux_dcop_fiesta_wpath.mat')
% % 18308, 18310 -> 37, 38
% % 63922,	63923/63941 -> 116, 117
% 
% load('/Volumes/TOSHIBA EXT/MINFLUX JS/241113/241113_05/241113-141106_minflux_dcop_555LP5p5_655_fbg70k_fiesta_wpath.mat')
% % 160, 161 -> 3, 4
% % 797, 798 -> 9, 10
% % 4391, 4393 -> 15, 16
% 
% load('/Volumes/TOSHIBA EXT/MINFLUX JS/241113/241113_05/241113-142248_minflux_555_bg30k_LP5p5_655_bg50k_LP8_fiesta_wpath.mat')
% % 11581, 11583 -> 28, 29
% 
% load('/Volumes/TOSHIBA EXT/MINFLUX JS/241113/241113_04/241113-125316_minflux_fiesta_wpath.mat')
% % 14860, 14862 -> 27, 28
% % 28708, 28714 -> 43, 44
% % 69431, 69432 -> 80, 81
% 
% load('/Volumes/TOSHIBA EXT/MINFLUX JS/241114/241114_03/241114-113643_minflux_fiesta_wpath.mat')
% % 67351, 67352 -> 307, 308
% 
% load('/Volumes/TOSHIBA EXT/MINFLUX JS/241115/241115_01/241115-084647_minflux_fiesta_wpath.mat')
% % 6313, 6314 -> 12, 13
% 
% load('/Volumes/TOSHIBA EXT/MINFLUX JS/241116/241116_15/241116-053824_minflux_555bg50k_655bg120k_fiesta_wpath.mat')
% % 42832, 42834 -> 251, 252




idx = [289, 290];
ch2correction = [0, 0];
Molecule1 = Molecule(idx(1)).Results;
Molecule2 = Molecule(idx(2)).Results;

% Define the range of time
time = linspace(0, 10, 100);  % Time points

% Define the x, y, and z (time) coordinates
x1 = Molecule1(:,3); x2 = Molecule2(:,3) + ch2correction(1);              % Example: x varies as sin(time)
y1 = Molecule1(:,4); y2 = Molecule2(:,4) + ch2correction(2);              % Example: y varies as cos(time)
z1 = Molecule1(:,2); z2 = Molecule2(:,2);                   % z-coordinate is time

% Plot the 3D line
figure;
plot3(x1, y1, z1, 'Color', Molecule(idx(1)).Color, 'LineWidth', 2);  % The 'LineWidth' argument makes the line thicker
hold on
plot3(x2, y2, z2, 'Color', Molecule(idx(2)).Color, 'LineWidth', 2);

% Add labels and title
xlabel('X');
ylabel('Y');
zlabel('Time');
title("3D Line Plot of " + Molecule(idx(1)).Name + " & " + Molecule(idx(2)).Name + " (X, Y, Time)");

% Adjust the view angle for better visualization
view(45, 30);
grid on;  % Add grid for better visualization


figure;
plot3(x1, y1, z1, 'Color', Molecule(idx(1)).Color, 'LineWidth', 2, 'LineStyle', '--');  % The 'LineWidth' argument makes the line thicker
scatter3(x1, y1, z1, 25, Molecule(idx(1)).Color, 'filled')
hold on
plot3(x2, y2, z2, 'Color', Molecule(idx(2)).Color, 'LineWidth', 2, 'LineStyle', '--');  % The 'LineWidth' argument makes the line thicker
scatter3(x2, y2, z2, 25, Molecule(idx(2)).Color, 'filled')

% Add labels and title
xlabel('X');
ylabel('Y');
zlabel('Time');
title("3D Line Plot of " + Molecule(idx(1)).Name + " & " + Molecule(idx(2)).Name + " (X, Y, Time)");

% Adjust the view angle for better visualization
view(45, 30);
grid on;  % Add grid for better visualization







%% ChatGPT

% % Define the range of time
% time = linspace(0, 10, 100);  % Time points
% 
% % Define the x, y, and z (time) coordinates
% x = sin(time);               % Example: x varies as sin(time)
% y = cos(time);               % Example: y varies as cos(time)
% z = time;                    % z-coordinate is time
% 
% % Plot the 3D line
% figure;
% plot3(x, y, z, 'LineWidth', 2);  % The 'LineWidth' argument makes the line thicker
% 
% % Add labels and title
% xlabel('X');
% ylabel('Y');
% zlabel('Time');
% title('3D Line Plot of (X, Y, Time)');
% 
% % Adjust the view angle for better visualization
% view(45, 30);
% grid on;  % Add grid for better visualization



%% Once we run path statistics, we should be able to do according to that path
% We have the information about raw xy and distance to the path. But we
% haven't drawn an identical path for both lines.

% But we could. We did something similar in the MAP9 handbook. Now I should
% do something similar.

% If we assume linear path, then this actually isn't too bad.

% Extrapolate linear path by beginning and end point, how far away in x and
% y we are for one plot
% Then calculate nearest neighbors to that line.

% ChatGPT again. It is 00:30 and I need to get up at 06:00.

% Define the points A, B, and P
% A = [1, 2];   % Point A on the line
A = Molecule(idx(1)).PathData(1,1:2);
B = Molecule(idx(1)).PathData(end,1:2);   % Point B on the line

% compute distance to line for idx(2) Molecule

xypath = zeros(length(Molecule(idx(2)).PathData(:,1)),2);
% spot for corrections if you have them from alignment
Molecule(idx(2)).PathData(:,1) = Molecule(idx(2)).PathData(:,1) + ch2correction(1);
Molecule(idx(2)).PathData(:,2) = Molecule(idx(2)).PathData(:,2) + ch2correction(2);

for j = 1:length(Molecule(idx(2)).PathData(:,1))
    P = Molecule(idx(2)).PathData(j,1:2); % External point P

    % Compute the vectors AB and AP
    AB = B - A;
    AP = P - A;
    
    % Compute the projection factor t
    t = dot(AP, AB) / dot(AB, AB);
    
    % Find the closest point C on the line
    % if t < 0
    %     C = A;  % Closest point is A
    % elseif t > 1
    %     C = B;  % Closest point is B
    % else
    C = A + t * AB;  % Point C is on the line segment
    % end
    
    % Compute the distance between P and C
    distance = norm(P - C);

    xypath(j,1:2) = [pdist2(A,C), distance];

end

figure
plot(Molecule(idx(1)).Results(:,2), Molecule(idx(1)).PathData(:,4), 'Color', Molecule(idx(1)).Color, 'LineWidth', 2)
hold on
plot(Molecule(idx(2)).Results(:,2), xypath(:,1), 'Color', Molecule(idx(2)).Color, 'LineWidth', 2)
title("Uncorrected on-axis path plot of " + Molecule(idx(1)).Name + " & " + Molecule(idx(2)).Name);

plot(Molecule(idx(1)).Results(:,2), Molecule(idx(1)).PathData(:,5), 'Color', [Molecule(idx(1)).Color, 0.2], 'LineWidth', 2)
plot(Molecule(idx(2)).Results(:,2), xypath(:,2), 'Color', [Molecule(idx(2)).Color, 0.2], 'LineWidth', 2)


% PRINT OUT INFO ABOUT TEMPORAL AND SPATIAL RESOLUTION
fprintf('Time and spatial resolution of ch1: \n')
mean(Molecule(idx(1)).Results(2:end,2) - Molecule(idx(1)).Results(1:end-1,2))
% mean(Molecule(idx(1)).Results(2:end,3) - Molecule(idx(1)).Results(1:end-1,3).^2 + Molecule(idx(1)).Results(2:end,4) - Molecule(idx(1)).Results(1:end-1,4).^2)

fprintf('Time and spatial resolution of ch2: \n')
mean(Molecule(idx(2)).Results(2:end,2) - Molecule(idx(2)).Results(1:end-1,2))
% mean(Molecule(idx(2)).Results(2:end,3) - Molecule(idx(2)).Results(1:end-1,3).^2 + Molecule(idx(2)).Results(2:end,4) - Molecule(idx(2)).Results(1:end-1,4).^2)

% P = [3, 3];   % External point P
% 
% % Compute the vectors AB and AP
% AB = B - A;
% AP = P - A;
% 
% % Compute the projection factor t
% t = dot(AP, AB) / dot(AB, AB);
% 
% % Find the closest point C on the line
% if t < 0
%     C = A;  % Closest point is A
% elseif t > 1
%     C = B;  % Closest point is B
% else
%     C = A + t * AB;  % Point C is on the line segment
% end
% 
% % Compute the distance between P and C
% distance = norm(P - C);
% 
% % Display the results
% fprintf('The closest point on the line is: [%.2f, %.2f, %.2f]\n', C);
% fprintf('The distance between the point and the line is: %.2f\n', distance);


% To do the actual after we have fit with FIONA
Ch1 = load('/Users/slivka/Documents/MATLAB/241109-151830_minflux_Joseph_interleaved_Cy3controller_bg50kCy3_100kCy5_fies/ 8842_fiona.mat');
Ch2 = load('/Users/slivka/Documents/MATLAB/241109-151830_minflux_Joseph_interleaved_Cy3controller_bg50kCy3_100kCy5_fies/ 8848_fiona.mat');

Ch1data = Ch1.data;
Ch2data = Ch2.data;

% figure()
nn1 = find(~isnan(Ch1data.xy(:,1)));
nn2 = find(~isnan(Ch2data.xy(:,1)));
% length(Ch1data.xy(nn,1));
figure()
plot(Ch1data.time, Ch1data.trace(nn1,3), 'k-', 'LineWidth', 1)
plot(Ch2data.time, Ch2data.trace(nn2,3)+88, 'r-', 'LineWidth', 1)



