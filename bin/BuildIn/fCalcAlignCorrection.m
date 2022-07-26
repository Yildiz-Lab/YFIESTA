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
for i = 1:length(ch1data.Molecule)
    mol = ch1data.Molecule(i);
        
    xbar = mean(mol.Results(:,3));
    ybar = mean(mol.Results(:,4));
    
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
                dir_store(k) = nan(1,2);
            end
        end
        
    end
        

end

%% calculate the mean separation of the molecules and plot a heatmap depending on their x,y positions
% Ideally, we would get a uniform heat map close to zero, but it is a real
% problem if the map is radial or has an asymmetric trend

mean_correction = mean(val_store,'omitnan') * mean(dir_store,'omitnan');

% do a 2d quiver plot to check that gradients are not too large (small
% spherical aberration


% make sure std is low enough that we can consider it even
std_mags = std(val_store.*dir_store, 'omitnan')

end

