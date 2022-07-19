function handles = ResetStepFitter(handles)
% RESETSTEPFITTER resets the parameters of the step-fitter, deletes the 
% current steps, and makes the Manual button invisible

% Written by Vladislav Belyy
% Last modified on 1/2/2012

% Erase previously fitted steps
handles.stepVector = 0;
set(handles.FitSteps, 'String', 'FIT');

set(handles.ManualStepFitting, 'Visible', 'off');
set(handles.ManualStepFitting, 'String', 'Adjust');

set(handles.SaveSteps, 'Visible', 'off');
set(handles.text72, 'Visible', 'off');
set(handles.AddDeleteStep, 'Visible', 'off');

iptPointerManager(handles.figureHandle, 'disable')




