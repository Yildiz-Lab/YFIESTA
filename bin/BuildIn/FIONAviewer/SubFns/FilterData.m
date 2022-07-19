function handles = FilterData(hObject, handles)
% Returns data filtered in accordance with given method and parameters
% All filtered data is saved in handles.PSD1Data_Long_Filt and
% handles.PSD1Data_Short_Filt
%
% Adapted from Tom Bilyard's original code by Vladislav Belyy


tic % start timing

WindowLength=str2double(get(handles.WindowLength,'string'));

FrameTime = str2double(get(handles.FrameLength, 'string'));
SR = 1000/FrameTime; % in Hz

x = RemoveNaNs(handles.PSD1Data_Long); % Long axis data
y = RemoveNaNs(handles.PSD1Data_Short); % Short axis data

Decimate=get(handles.Decimate,'Value');
Med=get(handles.MedianFilter,'Value');
Run=get(handles.RunningMean,'Value');
Butt=get(handles.Butterworth,'Value');
L1PWC=get(handles.L1PWC,'Value');

%         CurrentSelection='4: Filtered'
%         PrevSelection=handles.Ytype
%         PrevScale=handles.Yaxis

if Decimate 
    decimFactor = str2double(get(handles.DecimationFactor, 'String'));

    FilteredX_tmp = decimate(x, decimFactor);
    FilteredY_tmp = decimate(y, decimFactor);
    FilteredX = x;
    %keep length of decimated trace constant using linear interpolation
    for (i=2:length(FilteredX_tmp))
        for (j=1:decimFactor)
            FilteredX(i*decimFactor-j) = FilteredX_tmp(i) - j/decimFactor*(FilteredX_tmp(i)-FilteredX_tmp(i-1));
            FilteredY(i*decimFactor-j) = FilteredY_tmp(i) - j/decimFactor*(FilteredY_tmp(i)-FilteredY_tmp(i-1));
        end
    end
    FilteredX = FilteredX - (mean(FilteredX) - mean(x));
    FilteredY = FilteredY - (mean(FilteredY) - mean(y));
    
%     newSR = round(SR / decimFactor); % new sampling rate
% 
%     %Rebuild the time vector
%     handles.t_Filt = (1/newSR):(1/newSR): ...
%         (length(FilteredX)/newSR);

elseif Med %Median filtered
    FilteredX=medfilt1(x,WindowLength);
    FilteredY=medfilt1(y,WindowLength);
 
elseif Run % Running mean filter
    FilteredX=filter(ones(1,WindowLength)/(WindowLength),1,x);
    FilteredY=filter(ones(1,WindowLength)/(WindowLength),1,y);

elseif Butt % Butterworth filter
    CutoffFreq=str2double(get(handles.CutoffFreq,'string'))*2/SR;
    FilterOrder=str2double(get(handles.FilterOrder,'string'));
    [b,a] = butter(FilterOrder,CutoffFreq,'low');
    FilteredX=filter(b,a,x);
    FilteredY=filter(b,a,y);

elseif L1PWC % L1 piecewise-constant filter
    lambda=str2double(get(handles.lambda,'string'));
    FilteredX=l1tf_integ(x',lambda);
    FilteredY=l1tf_integ(y',lambda);

else
    FilteredX=x;
    FilteredY=y;

end

% save data
%FilteredX = FilteredX - mean(FilteredX) - mean(x);
%FilteredY = FilteredY - mean(FilteredY) - mean(x);

handles.PSD1Data_Long_Filt = FilteredX;
handles.PSD1Data_Short_Filt = FilteredY;

if ~Decimate
    handles.t_Filt = handles.t;
end

handles.FilteredFlag = 1;
filteringTime = toc; % end timing
disp(['Filtered in ', num2str(filteringTime), ' seconds']);

% Save the handles structure.
guidata(hObject,handles)