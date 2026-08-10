function app = AlignToMaxDx(app)

app = getappdata(app.fig, 'app');

window = app.FilterPanelFields.FilterWindow.Value;

if strcmp(app.FitPanelFields.DropFitAxis.Value,'Long-axis')
    offset = 0;
else
    offset = 1;
end
X = app.data.trace(:,[1,3,5]+offset);
x = X(:,1); xfit = X(:,2); cp = X(:,3);

% pad x with end values so that we don't have to worry about edge cases in taking mean differentials.
x = [x(1)*ones(2*window,1); x; x(end)*ones(2*window,1)];

% dx = abs(x(2:end) - x(1:end-1));

idx = find(cp)+2*window;
newcp = zeros(length(idx),1);

j=1;
for j = length(idx):-1:1
    i = idx(j);
    mn = i-window;
    mx = i+window;
    
    [mchange,midx,~] = largestMeanedChange(x(mn-window:mx+window),window); %similar to changepoint
    % [~, midx] = max(dx(mn:mx));
    % newcp(j) = i-window-1+midx;
    if mchange > 5
        newcp(j) = i-2*window-1+midx;
    else
        newcp(j) = [];
    end
    j=j-1;

end

newcp = newcp-2*window; %reset the index to be aligned with the original x
x = x(2*window+1:end-2*window); %Remove added x padding

% Now go through each new changepoint and refit the trace
cp = zeros(1,length(x)); cp(newcp) = 1;
newcp = [1; newcp; length(x)];

for i = 1:length(newcp)-1
    xfit(newcp(i):newcp(i+1)) = mean(x(newcp(i):newcp(i+1)),'omitnan');
end


% pass this into new stepFit changepoint
if strcmp(app.FitPanelFields.DropFitAxis.Value,'Long-axis')
    app.Data.stepVector = xfit'; % Update step vector
elseif strcmp(app.FitPanelFields.DropFitAxis.Value,'Short-axis')
    app.Data.shortstepVector = xfit'; % Update step vector
end
app = UserStepChangesTable(app, newcp, "align");
app = updateMainPlot(app);

stepFit.StepFit = xfit';
app = PackageTrace(app, stepFit);

setappdata(app.fig, 'app', app);



function [maxChange, idxMax, changes] = largestMeanedChange(x, halfWidth)
% largestMeanedChange
% Takes a 1D trace and computes the largest change between
% the mean of the left and right windows around each point
% (including the center point in both sides).
%
% INPUTS:
%   trace      - 1D vector
%   halfWidth  - default = 5 (for an 11-point window)
%
% OUTPUTS:
%   maxChange  - maximum absolute mean difference
%   idxMax     - index where this occurs
%   changes    - change value at every index

    x = x(:);   % ensure column vector
    n = length(x);

    changes = nan(n,1);

    for i = 1:n
        
        % Define window bounds safely
        leftStart  = max(1, i-halfWidth);
        rightEnd   = min(n, i+halfWidth);

        % Include center in both sides
        leftWindow  = x(leftStart:i);
        rightWindow = x(i:rightEnd);

        % Only compute if both exist
        if ~isempty(leftWindow) && ~isempty(rightWindow)
            leftMean  = mean(leftWindow);
            rightMean = mean(rightWindow);
            changes(i) = abs(rightMean - leftMean);
        end
    end

    [maxChange, idxMax] = max(changes, [], 'omitnan');
