function F = fCalcAlignCorrection(ch1file, ch2file)
%JS 2022/09/03
% Description:
%    function to correct mean offsets between two channels, ideally taken with
%    the same probes as to have 1 to 1 alignment
% Input:
%    ch1file Objects selected by uigetfile
%    ch2file Objects selected by uigetfile
% Output:
%    output will be x, y dependent translation correction for Ch2 to move
%    to Ch1 coordinates in nm

% too_large_nm = 1000; %FIONA using nanohole
too_large_nm = 20; %sqrt(3^2 + 7^2) %MINFLUX using 2C gold beads

% user choose option
if nargin < 1
% LoadDir = fShared('GetLoadDir');  % In FIESTA
LoadDir = uigetdir; % Outside of FIESTA (since we haven't quite gotten it in yet)
[baseName, folder] = uigetfile({'*.mat','FIESTA Data(*.mat)'},'Load FIESTA Objects for Channel 1',LoadDir,'MultiSelect','off');
ch1file = fullfile(folder, baseName);
[baseName, folder] = uigetfile({'*.mat','FIESTA Data(*.mat)'},'Load FIESTA Objects for Channel 2',LoadDir,'MultiSelect','off');
ch2file = fullfile(folder, baseName);
end

ch1data = load(ch1file);
ch2data = load(ch2file);

% JS Edit 2024/11/21 for MINFLUX data which already loads as molecules and
% has the channel corresponding
if ~isfield(ch1data,'Objects')

    % First a prefiltering/merging. If there are multiple Objects/Molecules within
    % a certain radius, we should combine them because they are actually the same localization.
    % This is because we don't want to mess up having pairings that look
    % like they are drifting away.
    
    xyc1 = nan(length(ch1data.Molecule),3);
    xyc2 = nan(length(ch2data.Molecule),3);
    for t = length(ch1data.Molecule):-1:1 %delete all ch1 files that don't belong to ch1
        if ch1data.Molecule(t).Channel == 1
            xyc1(t,1:2) = median(ch1data.Molecule(t).Results(:,3:4)); % package arrays for easy result comparison
            xyc1(t,3) = ch1data.Molecule(t).Channel;
        else
            xyc1(t,:) = [];
            ch1data.Molecule(t) = [];
        end
    end
    for t = length(ch2data.Molecule):-1:1 %delete all ch2 files that don't belong to ch2
        if ch2data.Molecule(t).Channel == 2
            xyc2(t,1:2) = median(ch2data.Molecule(t).Results(:,3:4)); % package arrays for easy result comparison
            xyc2(t,3) = ch2data.Molecule(t).Channel;
        else
            xyc2(t,:) = [];
            ch2data.Molecule(t) = [];
        end
    end

    R = 150; %nm

    for n = length(ch1data.Molecule):-1:1
        mol_combine = find(pdist2(xyc1(n,1:2), xyc1(:,1:2)) < R);
        mol_combine(end) = []; %erase n match, or else it will delete itself
        if ~isempty(mol_combine)
            ch1data.Molecule(mol_combine(end)).Results = [ch1data.Molecule(mol_combine(end)).Results(:,:); ch1data.Molecule(n).Results(:,:)];
            xyc1(mol_combine(end),1:2) = median(ch1data.Molecule(mol_combine(end)).Results(3:4));
            ch1data.Molecule(n) = []; % Delete
            xyc1(n,:) = [];
        end
    end
    for n = length(ch2data.Molecule):-1:1
        mol_combine = find(pdist2(xyc2(n,1:2), xyc2(:,1:2)) < R);
        mol_combine(end) = []; %erase n match, or else it will delete itself
        if ~isempty(mol_combine)
            ch2data.Molecule(mol_combine(end)).Results = [ch2data.Molecule(mol_combine(end)).Results(:,:); ch2data.Molecule(n).Results(:,:)];
            xyc2(mol_combine(end),1:2) = median(ch2data.Molecule(mol_combine(end)).Results(3:4));
            ch2data.Molecule(n) = []; % Delete
            xyc2(n,:) = [];
        end
    end
    fprintf(['xyc', num2str(length(xyc1)), '\n'])
    fprintf(['xyc', num2str(length(xyc2)), '\n'])
    

    ch1data.Objects = {};
    ch2data.Objects = {};
    
    ch1cnt = 0;
    ch2cnt = 0;

    for s = 1:length(ch1data.Molecule)
        ch1cnt = ch1cnt + 1;
        % ch1data.Objects{ch1cnt} = struct();
        ch1data.Objects{1}.center_x(ch1cnt) = mean(ch1data.Molecule(s).Results(:,3));
        ch1data.Objects{1}.center_y(ch1cnt) = mean(ch1data.Molecule(s).Results(:,4));
    end


    for s = 1:length(ch2data.Molecule)
        ch2cnt = ch2cnt + 1;
        % ch2data.Objects{1} = struct();
        ch2data.Objects{1}.center_x(ch2cnt) = mean(ch2data.Molecule(s).Results(:,3));
        ch2data.Objects{1}.center_y(ch2cnt) = mean(ch2data.Molecule(s).Results(:,4));
    end

    % for s = 1:length(ch1data.Molecule)
    %     if ch1data.Molecule(s).Channel == 1
    %         ch1cnt = ch1cnt + 1;
    %         % ch1data.Objects{ch1cnt} = struct();
    %         ch1data.Objects{1}.center_x(ch1cnt) = median(ch1data.Molecule(s).Results(:,3));
    %         ch1data.Objects{1}.center_y(ch1cnt) = median(ch1data.Molecule(s).Results(:,4));
    %     elseif ch1data.Molecule(s).Channel == 2
    %         ch2cnt = ch2cnt + 1;
    %         % ch2data.Objects{ch2cnt} = struct();
    %         ch2data.Objects{1}.center_x(ch2cnt) = median(ch2data.Molecule(s).Results(:,3));
    %         ch2data.Objects{1}.center_y(ch2cnt) = median(ch2data.Molecule(s).Results(:,4));
    %     end
    % 
    % 
    % end
    
    % % if they are actually different files, so then we should separate.
    % if ch1file ~= ch2file
    %     ch2data.Objects = {};
    %     ch2cnt = 0;
    %     for s = 1:length(ch2data.Molecule)
    %         ch2cnt = ch2cnt + 1;
    %         ch2data.Objects{ch2cnt} = struct();
    %         ch2data.Objects{ch2cnt}.center_x = median(ch2data.Molecule(s).Results(:,3));
    %         ch2data.Objects{ch2cnt}.center_y = median(ch2data.Molecule(s).Results(:,4));
    %     end
    % end

    
end


%% Make an array to store all matches from all frames
% We are going to make a map that takes Ch2 into Ch1
% [Ch1 center_x, Ch1 center_y, Ch2 center_x, Ch2 center_y, xcorrect, ycorrect]

GridPoints = [];
m=0;

%% Calculate each frame of objects separately, but will stitch together before end

for k = 1:length(ch2data.Objects)

    frobj1 = ch1data.Objects{k};
    frobj2 = ch2data.Objects{k};
    
    %frobj1.center_x
    
    %% Get centers for Ch2 for time save
    Ch2 = zeros(length(frobj2.center_x),2);
    for j = 1:length(frobj2.center_x)
        Ch2(j,1) = frobj2.center_x(j);
        Ch2(j,2) = frobj2.center_y(j);
    end

    %% Go through and try to pair each molecule with another
    idx_store = nan(length(frobj1.center_x),1);
    val_store = nan(length(frobj1.center_x),1);
    dir_store = nan(length(frobj1.center_x),2);
    Ch1 = zeros(length(frobj1.center_x),2);
    for i = 1:length(frobj1.center_x)
        
        xbar = frobj1.center_x(i);
        ybar = frobj1.center_y(i);
        Ch1(i,:) = [xbar, ybar];

        % magnitude
        [val, idx] = min(sqrt((Ch2(:,1)-xbar).^2 + (Ch2(:,2)-ybar).^2));

        if val < too_large_nm
            % make sure this molecule hasn't already been paired
            if ~ismember(idx_store, idx)
                idx_store(i) = idx;
                val_store(i) = val;
                dir_store(i,:) = [Ch2(idx,1)-xbar,Ch2(idx,2)-ybar]/val;
            % if it has, then is the new molecule closer. Else just keep it nan
            else
                [~,k] = find(idx_store == idx);
                if val < val_store(k)
                    idx_store(i) = idx;
                    val_store(i) = val;
                    dir_store(i,:) = [Ch2(idx,1)-xbar,Ch2(idx,2)-ybar]/val;
                    idx_store(k) = NaN;
                    val_store(k) = NaN;
                    dir_store(k,:) = nan(1,2);
                end
            end

        end


    end
    
    %% we have now paired all frame molecules. We should append this to a running list
    
    %Just in case one wants to see the percentage that didn't get matched
    % If this number is too high, valuable x,y information may be missing
    m = m + length(idx_store); 
    mnan = find(~isnan(idx_store));
    
    GridPoints = [GridPoints; Ch1(mnan,1), Ch1(mnan,2), Ch2(idx_store(mnan),1), Ch2(idx_store(mnan),2), val_store(mnan).*dir_store(mnan,1), val_store(mnan).*dir_store(mnan,2)];
    
    
end


fprintf("Total Grid Points Matched: " + num2str(size(GridPoints,1)) + "\n")
fprintf("Found Objects Matched: " + num2str(round(size(GridPoints,1)*100/m,2)) + " percent")
fprintf("\n Domain of Matching: x [" + min(GridPoints(3,:)) + "," + max(GridPoints(3,:)) + "]  y [" + min(GridPoints(4,:)) + "," + max(GridPoints(4,:)) + "] \n")

sum(isnan(GridPoints(:,5)))

fprintf("Before Interpolation Correction Error: x: " + num2str(round(mean(GridPoints(:,5)),2)) + "  y: " + num2str(round(mean(GridPoints(:,6)),2)) + "\n")

deltax = fCubicSplineInterp(GridPoints(:,3), GridPoints(:,4), GridPoints(:,5));
deltay = fCubicSplineInterp(GridPoints(:,3), GridPoints(:,4), GridPoints(:,6));
fprintf("After Interpolation Correction Error: x: " + num2str(round(deltax,2)) + "  y: " + num2str(round(deltay,2)) + "\n")

%% Now that it is tested, get the function and eventually export it to the mainGui

dx = fCubicSplineInterp(GridPoints(:,3), GridPoints(:,4), GridPoints(:,5),1);
dy = fCubicSplineInterp(GridPoints(:,3), GridPoints(:,4), GridPoints(:,6),1);
F = {dx, dy};

% So now to get a new point, you use
% x2 = 30000; y2 = 15000;
% xy2prime = [x2 + dx(x2,y2), y2 + dy(x2,y2)]

%% OLD CODE %%

% %% Calculate the means for Ch2 for time save
% Ch2 = zeros(length(ch2data.Molecule),2);
% for j = 1:length(ch2data.Molecule)
%     for mol = ch2data.Molecule(j)
%         Ch2(j,1) = mean(mol.Results(:,3));
%         Ch2(j,2) = mean(mol.Results(:,4));
%     end
% end
% 
% %% Go through and try to pair each molecule with another
% idx_store = nan(length(ch1data.Molecule),1);
% val_store = nan(length(ch1data.Molecule),1);
% dir_store = nan(length(ch1data.Molecule),2);
% Ch1 = zeros(length(ch1data.Molecule),2);
% for i = 1:length(ch1data.Molecule)
%     mol = ch1data.Molecule(i);
%         
%     xbar = mean(mol.Results(:,3));
%     ybar = mean(mol.Results(:,4));
%     Ch1(i,:) = [xbar, ybar];
%     
%     % magnitude
%     [val, idx] = min(sqrt((Ch2(:,1)-xbar).^2 + (Ch2(:,2)-ybar).^2));
% 
%     if val < too_large_nm
%         % make sure this molecule hasn't already been paired
%         if ~ismember(idx_store, idx)
%             idx_store(i) = idx;
%             val_store(i) = val;
%             dir_store(i,:) = [Ch2(idx,1)-xbar,Ch2(idx,2)-ybar]/val;
%         % if it has, then is the new molecule closer. Else just keep it nan
%         else
%             [~,k] = find(idx_store == idx);
%             if val < val_store(k)
%                 idx_store(i) = idx;
%                 val_store(i) = val;
%                 dir_store(i,:) = [Ch2(idx,1)-xbar,Ch2(idx,2)-ybar]/val;
%                 idx_store(k) = NaN;
%                 val_store(k) = NaN;
%                 dir_store(k,:) = nan(1,2);
%             end
%         end
%         
%     end
%         
% 
% end

% mask out nan
% mnan = ~isnan(dir_store(:,1));
% 
% z = val_store(mnan,1).*dir_store(mnan,:);
% 
% deltax = fCubicSplineInterp(Ch2(idx_store(mnan),1), Ch2(idx_store(mnan),2), z(:,1))
% deltay = fCubicSplineInterp(Ch2(idx_store(mnan),1), Ch2(idx_store(mnan),2), z(:,2))

%% calculate the mean separation of the molecules and plot a heatmap depending on their x,y positions
% Ideally, we would get a uniform heat map close to zero, but it is a real
% problem if the map is radial or has an asymmetric trend


% % get an initial guess for the center from LSS
% A = [ -dir_store(mnan,2)./dir_store(mnan,1), ones(sum(mnan),1)];
% B = Ch2(idx_store(mnan),2) - dir_store(mnan,2) ./ dir_store(mnan,1) .* Ch2(idx_store(mnan),1);
% LSScenter = A\B;
% 
% % https://www.mathworks.com/matlabcentral/answers/119001-2d-data-fitting-surface
% radfit = @(C,XY) C(1) .* sqrt((XY(:,1)-C(2)).^2 + (XY(:,2)-C(3)).^2);
% C = lsqcurvefit(radfit, [1 LSScenter(1) LSScenter(2)], Ch2(idx_store(mnan),:), val_store(mnan));
% 
% % res = 1000
% % want to fit on a meshgrid
% % X = 0:res:88722;
% % Y = 0:res:43200;
% % 
% % XY = zeros(length(X)*length(Y),2);
% % for i = 1:length(X)
% %     for j = 1:length(Y)
% %         XY((i-1)*length(Y)+j,:) = [X(i), Y(j)];
% %     end
% % end
% 
% %reshape(radfit(C,XY),length(Y),length(X))
% 
% figure(88)
% hold on
% scatter(Ch1(:,1),Ch1(:,2),'r','DisplayName','Ch1')
% scatter(Ch2(:,1),Ch2(:,2),'g','DisplayName','Ch2')
% uv = [Ch2(:,1)-C(2), Ch2(:,2)-C(3)];
% uv = radfit(C,Ch2).*uv./vecnorm(uv')';
% 
% % add in the direction since we only deal with magnitude and the code
% % doesn't know which way to apply
% 
% pd = mean(vecnorm(Ch1(mnan,:)' - Ch2(idx_store(mnan),:)'+uv(idx_store(mnan),:)'));
% md = mean(vecnorm(Ch1(mnan,:)' - Ch2(idx_store(mnan),:)'-uv(idx_store(mnan),:)'));
% if pd > md
%     uv = -uv;
% end
% 
% scatter(Ch2(:,1)-uv(:,1), Ch2(:,2)-uv(:,2),'b','DisplayName','Ch2 Corrected')
% legend()
% 
% scatter(C(2), C(3), 100, 'filled')
% 
% fprintf(strcat('Mean Distance Off:  ', num2str(mean(vecnorm(Ch1(mnan,:)' - Ch2(idx_store(mnan),:)'))),'\n'))
% fprintf(strcat('Mean Distance Off following radial correction:  ', num2str(mean(vecnorm(Ch1(mnan,:)' - Ch2(idx_store(mnan),:)'+uv(idx_store(mnan),:)'))),'\n'))
% 
% 
% mean_correction = mean(val_store,'omitnan') * mean(dir_store,'omitnan');
% scatter(Ch2(:,1)-mean_correction(1), Ch2(:,2)-mean_correction(2),'k','DisplayName','Translation Ch2 Corrected')
% OOON = ones(length(idx_store(mnan)),1); % I was sick of writing ones(length(idx_store(mnan)),1)
% fprintf(strcat('Mean Distance Off following translational correction:  ', num2str(mean(vecnorm(Ch1(mnan,:)' - (Ch2(idx_store(mnan),:) - [mean_correction(1)*OOON, mean_correction(2)*OOON])' ))),'\n'))
% 
% % do a 2d quiver plot to check that gradients are not too large (small
% % spherical aberration
% 
% 
% % make sure std is low enough that we can consider it even
% std_mags = std(val_store.*dir_store, 'omitnan');

end

