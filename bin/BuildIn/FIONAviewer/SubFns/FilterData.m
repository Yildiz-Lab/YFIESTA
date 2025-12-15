function app = FilterData(app)
% Returns data filtered in accordance with given method and parameters
% All filtered data is saved in handles.PSD1Data_Long_Filt and
% handles.PSD1Data_Short_Filt
%
% Adapted from Tom Bilyard's original code by Vladislav Belyy


WindowLength=app.FilterPanelFields.FilterWindow.Value;

% x = app.Data.t; % Long axis data
% if strcmp(app.FitPanelFields.DropFitAxis.Value,'Long-axis')
%     y = app.Data.PSD1Data_Long;
% else
%     y = app.Data.PSD1Data_Short;
% end
x = app.Data.PSD1Data_Long;
y = app.Data.PSD1Data_Short;

% app.FilterPanelFields.FilterGroup.SelectedObject.Text
if strcmp(app.FilterPanelFields.FilterGroup.SelectedObject.Text,'Decimate')
    t = 1:length(x);
    t_ds = 1:WindowLength:length(x);
    FilteredX = interp1(t, x, t_ds, 'linear');
    FilteredY = interp1(t, y, t_ds, 'linear');

elseif strcmp(app.FilterPanelFields.FilterGroup.SelectedObject.Text,'Median Filter') %Median filtered
    FilteredX=medfilt1(x,WindowLength);
    FilteredY=medfilt1(y,WindowLength);
 
elseif strcmp(app.FilterPanelFields.FilterGroup.SelectedObject.Text,'Running Mean') % Running mean filter
    FilteredX=filter(ones(1,WindowLength)/(WindowLength),1,x);
    FilteredY=filter(ones(1,WindowLength)/(WindowLength),1,y);

elseif strcmp(app.FilterPanelFields.FilterGroup.SelectedObject.Text,'Butterworth') % Butterworth filter
    % CutoffFreq=str2double(get(handles.CutoffFreq,'string'))*2/SR;
    % FilterOrder=str2double(get(handles.FilterOrder,'string'));
    % [b,a] = butter(FilterOrder,CutoffFreq,'low');
    % FilteredX=filter(b,a,x);
    % FilteredY=filter(b,a,y);
    order = app.FilterPanelFields.FilterOrder.Value;  % filter order (4 is usually fine)
    Wn = 1 / WindowLength;  % normalized cutoff (1/downsample rate) so it is below new Nyquist
    
    [b,a] = butter(order, Wn); % design lowpass
    x_filt = filtfilt(b,a,x);  % zero-phase filtering
    y_filt = filtfilt(b,a,y);
    
    % subsample
    t_ds = 1:WindowLength:length(x_filt);
    FilteredX = x_filt(t_ds);
    FilteredY = y_filt(t_ds);
    


elseif strcmp(app.FilterPanelFields.FilterGroup.SelectedObject.Text,'L1 Piecewise Constant') % L1 piecewise-constant filter
    % xmin​∥y−x∥22​+λ∥Dx∥1
    % Higher lambda penalizes more steps. i.e. high lambda -> heavy
    % smoothing -> fewer steps
    lambda=app.FilterPanelFields.FilterLambdaPenalty.Value;
    % FilteredX=l1tf_integ(x',lambda);
    % FilteredY=l1tf_integ(y',lambda);
    FilteredX = l1tf_stepfilter(x', lambda);
    FilteredY = l1tf_stepfilter(y', lambda);

else
    FilteredX=x;
    FilteredY=y;
end

app.Data.PSD1Data_Long_Filt = FilteredX;
app.Data.PSD1Data_Short_Filt = FilteredY;

if strcmp(app.FilterPanelFields.FilterGroup.SelectedObject.Text,'Decimate')
    app.Data.t_Filt = app.Data.t(t_ds);
elseif strcmp(app.FilterPanelFields.FilterGroup.SelectedObject.Text,'Butterworth')
    app.Data.t_Filt = app.Data.t(t_ds);
else
    app.Data.t_Filt = app.Data.t;
end

app.Data.FilteredFlag = 1;



function xhat = l1tf_stepfilter(x, lambda) %chatGPT
%L1TF_STEPFILTER  Piecewise-constant L1 trend filtering (TVD)
%
%   xhat = l1tf_stepfilter(x, lambda)
%
%   Inputs:
%       x      - input signal (row or column vector)
%       lambda - regularization strength
%
%   Output:
%       xhat   - piecewise-constant filtered signal
%
%   This implements fast 1D Total Variation Denoising
%   using Condat's O(N) algorithm:
%       xhat = argmin_z  0.5*||z - x||^2 + lambda * TV(z)
%
%   TV(z) = sum |z(i+1) - z(i)|, promoting flat regions (steps).
%
%   A larger lambda = fewer, bigger steps.
%   A smaller lambda = keeps more detail.
%
%   No toolboxes required.

x = x(:);                 % ensure column vector
N = length(x);

% Initialize
xhat = zeros(N,1);
k = 1; k0 = 1;
vmin = x(1) - lambda;
vmax = x(1) + lambda;
umin = lambda - x(1);
umax = -lambda - x(1);

for i = 2:N
    % Update bounds
    vmin = min(vmin, x(i)-lambda);
    vmax = max(vmax, x(i)+lambda);
    umin = umin + x(i) - lambda;
    umax = umax + x(i) + lambda;

    % Check feasibility
    if vmin > umax
        for j = k0:k
            xhat(j) = vmin;
        end
        k0 = k+1;
        k = k0;
        vmin = x(i)-lambda;
        vmax = x(i)+lambda;
        umin = lambda-x(i);
        umax = -lambda-x(i);

    elseif vmax < umin
        for j = k0:k
            xhat(j) = vmax;
        end
        k0 = k+1;
        k = k0;
        vmin = x(i)-lambda;
        vmax = x(i)+lambda;
        umin = lambda-x(i);
        umax = -lambda-x(i);

    else
        k = k + 1;
    end
end

% Final averaging step
v = 0.5*(vmin + vmax);
for j = k0:N
    xhat(j) = v;
end




