function mean_correction = fCalcAlignCorrection
%JS 2022/07/26
% Description:
%    function to correct mean offsets between two channels, ideally taken with
%    the same probes as to have 1 to 1 alignment
% Input:
%    ch1file Molecules selected by uigetfile
%    ch2file Molecules selected by uigetfile
% Output:
%    output will be an x, y translation correction for Ch2 to move on top of
%    Ch1 in nm

too_large_nm = 3000;

% user choose option
% LoadDir = fShared('GetLoadDir');  % In FIESTA
LoadDir = uigetdir; % Outside of FIESTA (since we haven't quite gotten it in yet)
[baseName, folder] = uigetfile({'*.mat','FIESTA Data(*.mat)'},'Load FIESTA Tracks Channel 1',LoadDir,'MultiSelect','off');
ch1file = fullfile(folder, baseName);
[baseName, folder] = uigetfile({'*.mat','FIESTA Data(*.mat)'},'Load FIESTA Tracks Channel 2',LoadDir,'MultiSelect','off');
ch2file = fullfile(folder, baseName);

% ch1file = '.mat';
% ch2file = '.mat';

ch1data = load(ch1file);
ch2data = load(ch2file);

%% Calculate the means for Ch2 for time save
Ch2 = zeros(length(ch2data.Molecule),2);
for j = 1:length(ch2data.Molecule)
    for mol = ch2data.Molecule(j)
        Ch2(j,1) = mean(mol.Results(:,3));
        Ch2(j,2) = mean(mol.Results(:,4));
    end
end

%% Go through and try to pair each molecule with another
idx_store = nan(length(ch1data.Molecule),1);
val_store = nan(length(ch1data.Molecule),1);
dir_store = nan(length(ch1data.Molecule),2);
Ch1 = zeros(length(ch1data.Molecule),2);
for i = 1:length(ch1data.Molecule)
    mol = ch1data.Molecule(i);
        
    xbar = mean(mol.Results(:,3));
    ybar = mean(mol.Results(:,4));
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

%% calculate the mean separation of the molecules and plot a heatmap depending on their x,y positions
% Ideally, we would get a uniform heat map close to zero, but it is a real
% problem if the map is radial or has an asymmetric trend



% mask out nan
mnan = ~isnan(dir_store(:,1));
% get an initial guess for the center from LSS
A = [ -dir_store(mnan,2)./dir_store(mnan,1), ones(sum(mnan),1)];
B = Ch2(idx_store(mnan),2) - dir_store(mnan,2) ./ dir_store(mnan,1) .* Ch2(idx_store(mnan),1);
LSScenter = A\B;

% https://www.mathworks.com/matlabcentral/answers/119001-2d-data-fitting-surface
radfit = @(C,XY) C(1) .* sqrt((XY(:,1)-C(2)).^2 + (XY(:,2)-C(3)).^2);
C = lsqcurvefit(radfit, [1 LSScenter(1) LSScenter(2)], Ch2(idx_store(mnan),:), val_store(mnan));

% res = 1000
% want to fit on a meshgrid
% X = 0:res:88722;
% Y = 0:res:43200;
% 
% XY = zeros(length(X)*length(Y),2);
% for i = 1:length(X)
%     for j = 1:length(Y)
%         XY((i-1)*length(Y)+j,:) = [X(i), Y(j)];
%     end
% end

%reshape(radfit(C,XY),length(Y),length(X))

figure(88)
hold on
scatter(Ch1(:,1),Ch1(:,2),'r','DisplayName','Ch1')
scatter(Ch2(:,1),Ch2(:,2),'g','DisplayName','Ch2')
uv = [Ch2(:,1)-C(2), Ch2(:,2)-C(3)];
uv = radfit(C,Ch2).*uv./vecnorm(uv')';
scatter(Ch2(:,1)+uv(:,1), Ch2(:,2)+uv(:,2),'b','DisplayName','Ch2 Corrected')
legend()

fprintf(strcat('Mean Distance Off:  ', num2str(mean(vecnorm(Ch1(mnan,:)' - Ch2(idx_store(mnan),:)'))),'\n'))
fprintf(strcat('Mean Distance Off following radial correction:  ', num2str(mean(vecnorm(Ch1(mnan,:)' - Ch2(idx_store(mnan),:)'-uv(idx_store(mnan),:)'))),'\n'))


mean_correction = mean(val_store,'omitnan') * mean(dir_store,'omitnan');
scatter(Ch2(:,1)-mean_correction(1), Ch2(:,2)+uv(:,2)-mean_correction(2),'k','DisplayName','Translation Ch2 Corrected')

% do a 2d quiver plot to check that gradients are not too large (small
% spherical aberration


% make sure std is low enough that we can consider it even
std_mags = std(val_store.*dir_store, 'omitnan');

end

