function app = UserStepChangesTable(app, stepIdx, mode, varargin)
% UserStepChangesTable  Create/append/align step-change rows, removing duplicates.
% Modified to preserve existing Action and String when aligning.

% Validate and normalize step indices
if isempty(stepIdx)
    return
end
stepIdx = double(stepIdx(:)); % column vector

% Parse name-value inputs
p = inputParser;
addParameter(p, "Action", "", @(x) ischar(x) || isstring(x));
addParameter(p, "UserString", "", @(x) ischar(x) || isstring(x));
parse(p, varargin{:});
actionVal = string(p.Results.Action);
userStr   = string(p.Results.UserString);

% Determine Axis from dropdown (safely)
axisVal = "unknown";
try
    axisVal = string(app.FitPanelFields.DropFitAxis.Value);
catch
end

% Determine Action and String based on mode
switch string(mode)
    case "fit"
        actionVal = "add";
        userStr = "fit";
    case "user"
        if actionVal == ""
            % try
            %     if app.FitPanelFields.BtnDelete.UserData.pressed
            %         actionVal = "delete";
            %     end
            % catch
            % end
            try
                if app.FitPanelFields.BtnAdd.UserData.pressed
                    actionVal = "add";
                end
            catch
            end
            if actionVal == ""
                actionVal = "add";
            end
        end
        if userStr == ""
            try
                userStr = string(app.FitPanelFields.DropFitString.Value);
            catch
                userStr = "";
            end
        end
    case "align"
        if actionVal == ""
            % try
            %     if app.FitPanelFields.BtnDelete.UserData.pressed
            %         actionVal = "delete";
            %     end
            % catch
            % end
            try
                if app.FitPanelFields.BtnAdd.UserData.pressed
                    actionVal = "add";
                end
            catch
            end
            if actionVal == ""
                actionVal = "add";
            end
        end
        if userStr == ""
            try
                userStr = string(app.FitPanelFields.DropFitString.Value);
            catch
                userStr = "";
            end
        end
    otherwise
        error("Mode must be 'fit', 'user', or 'align'.")
end

% Determine filter window (fallback to 0 if missing or invalid)
try
    fw = app.FilterPanelFields.FilterWindow.Value;
    filterWindow = double(fw);
    if ~isfinite(filterWindow) || filterWindow < 0
        filterWindow = 0;
    end
catch
    filterWindow = 0;
end

% Initialize StepChangesTable if empty / missing
if ~isfield(app.Data, "StepChangesTable") || isempty(app.Data.StepChangesTable)
    T = table( ...
        double.empty(0,1), string.empty(0,1), string.empty(0,1), string.empty(0,1), ...
        'VariableNames', {'Index','Axis','Action','String'} );
else
    T = app.Data.StepChangesTable;
    if ~istable(T)
        try
            T = struct2table(T);
        catch
            try
                T = cell2table(T, 'VariableNames', {'Index','Axis','Action','String'});
            catch
                T = table([],[],[],[], 'VariableNames', {'Index','Axis','Action','String'});
            end
        end
    end
    missing = setdiff({'Index','Axis','Action','String'}, T.Properties.VariableNames);
    for k = 1:numel(missing)
        T.(missing{k}) = repmat("", height(T), 1);
    end
end

% Normalize columns
if ~isnumeric(T.Index)
    try
        T.Index = double(T.Index);
    catch
        T.Index = nan(height(T),1);
    end
end
T.Axis = string(T.Axis);
T.Action = string(T.Action);
T.String = string(T.String);

switch string(mode)
    case "user"
        % Remove existing rows within filterWindow of any new stepIdx (same Axis)
        if ~isempty(T) && ~isempty(stepIdx)
            existingIdx = T.Index(:);
            existingAxis = T.Axis(:);
            keepMask = true(height(T),1);
            if filterWindow > 0
                for i = 1:height(T)
                    if existingAxis(i) == axisVal
                        if any(abs(existingIdx(i) - stepIdx) <= filterWindow)
                            keepMask(i) = false;
                        end
                    end
                end
            else
                sameAxisMask = (existingAxis == axisVal);
                dupMask = ismember(existingIdx, stepIdx) & sameAxisMask;
                keepMask(dupMask) = false;
            end
            T = T(keepMask, :);
        end
        nNew = numel(stepIdx);
        Tnew = table(stepIdx, repmat(axisVal,nNew,1), repmat(actionVal,nNew,1), repmat(userStr,nNew,1), ...
            'VariableNames', {'Index','Axis','Action','String'});
        T = [T; Tnew];
    case "fit"
        nNew = numel(stepIdx);
        Tnew = table(stepIdx, repmat(axisVal,nNew,1), repmat(actionVal,nNew,1), repmat(userStr,nNew,1), ...
            'VariableNames', {'Index','Axis','Action','String'});
        if ~isempty(T)
            existingIdx = T.Index(:);
            existingAxis = T.Axis(:);
            keepMask = true(height(T),1);
            if filterWindow > 0
                for i = 1:height(T)
                    if existingAxis(i) == axisVal
                        if any(abs(existingIdx(i) - stepIdx) <= filterWindow)
                            keepMask(i) = false;
                        end
                    end
                end
            else
                sameAxisMask = (existingAxis == axisVal);
                dupMask = ismember(existingIdx, stepIdx) & sameAxisMask;
                keepMask(dupMask) = false;
            end
            T = T(keepMask,:);
        end
        T = [T; Tnew];
    case "align"
        % Preserve existing Action and String when aligning.
        if isempty(T)
            % if actionVal ~= "delete"
            %     nNew = numel(stepIdx);
            %     Tnew = table(stepIdx, repmat(axisVal,nNew,1), repmat(actionVal,nNew,1), repmat(userStr,nNew,1), ...
            %         'VariableNames', {'Index','Axis','Action','String'});
            %     T = [T; Tnew];
            % end
        else
            for k = 1:numel(stepIdx)
                newIdx = stepIdx(k);
                matchMask = (T.Axis == axisVal) & (abs(T.Index - newIdx) <= filterWindow);
                matches = find(matchMask);
                if isempty(matches)
                    % if actionVal ~= "delete"
                    %     % append new row using provided action and string
                    %     T = [T; table(newIdx, axisVal, actionVal, userStr, 'VariableNames', T.Properties.VariableNames)];
                    % end
                else
                    % if actionVal == "delete"
                    %     % remove matched rows
                    %     T(matches,:) = [];
                    % else
                        % update matched rows' Index and Axis only; preserve Action and String
                        for m = 1:numel(matches)
                            row = matches(m);
                            T.Index(row) = newIdx;
                            T.Axis(row) = axisVal;
                            % Do NOT overwrite T.Action(row) or T.String(row)
                        end
                    % end
                end
            end
            % Remove exact duplicates that may have been created (keep first occurrence)
            if height(T) > 1
                key = strcat(string(T.Index), "|", T.Axis, "|", T.Action, "|", T.String);
                [~, ia] = unique(key, 'stable');
                maskKeep = false(height(T),1);
                maskKeep(ia) = true;
                T = T(maskKeep, :);
            end
        end
    otherwise
        error("Unknown mode.")
end

% Final sort by Index
try
    T = sortrows(T, 'Index');
catch
end

% Write back into app.Data.StepChangesTable and app.StepChangesTable (if present)
app.Data.StepChangesTable = T;
try
    app.StepChangesTable = T;
catch
end

T

% Persist app handle if original code expects it
try
    setappdata(app.fig, 'app', app);
catch
end

end






% function app = UserStepChangesTable(app, stepIdx, mode, varargin)
% % UserStepChangesTable  Create/append step-change rows, removing duplicates.
% %
% % UserStepChangesTable(app, stepIdx, mode)
% %   mode = "fit"  -> Action="add", String="fit"
% %   mode = "user" -> Action and UserString should be provided via name-value:
% %       UpdateStepChangesTable(..., "Action", actionValue, "UserString", userStr)
% %   If Action is not provided in "user" mode, the function will attempt to
% %   infer it from toggle/button Value properties:
% %       app.FitPanelFields.BtnAdd.Value or app.FitPanelFields.BtnDelete.Value
% %
% % Inputs:
% %   app       - app object
% %   stepIdx   - numeric scalar or vector of step indices to add/append
% %   mode      - "fit" or "user"
% %   Name-Value:
% %       "Action"     - "add" or "delete" (string/char)
% %       "UserString" - string to place in the String column (default: "")
% %
% % The function updates:
% %   app.FilterPanelFields.FilterWindow         
% %   app.StepChangesTable                      (table copy)
% 
% % Validate and normalize step indices
% if isempty(stepIdx)
%     return
% end
% stepIdx = double(stepIdx(:)); % column vector
% 
% % Parse name-value inputs
% p = inputParser;
% addParameter(p, "Action", "", @(x) ischar(x) || isstring(x));
% addParameter(p, "UserString", "", @(x) ischar(x) || isstring(x));
% parse(p, varargin{:});
% actionVal = string(p.Results.Action);
% userStr   = string(p.Results.UserString);
% 
% % Determine Axis from dropdown (safely)
% try
%     axisVal = string(app.FitPanelFields.DropFitAxis.Value);
% catch
%     axisVal = "unknown";
% end
% 
% % Determine Action and String based on mode
% switch string(mode)
%     case "fit"
%         actionVal = "add";
%         userStr = "fit";
%     case "user"
%         % If Action not provided, try to infer from BtnAdd/BtnDelete Value
%         if actionVal == ""
% 
%             if app.FitPanelFields.BtnDelete.UserData.pressed
%                 actionVal = "delete";
%             end
% 
%             if app.FitPanelFields.BtnAdd.UserData.pressed
%                 actionVal = "add";
%             end
% 
%             % fallback
%             if actionVal == ""
%                 actionVal = "add";
%             end
%         end
%         if userStr == ""
%             % Try to read user edit field if it exists
% 
%             try
%                 userStr = string(app.FitPanelFields.DropFitString.Value);
%             catch
%                 userStr = "";
%             end
% 
%         end
%     otherwise
%         error("Mode must be 'fit' or 'user'.")
% end
% 
% 
% 
% % Determine filter window (fallback to 0 if missing or invalid)
% try
%     filterWindow = double(app.FilterPanelFields.FilterWindow.Value);
%     if ~isfinite(filterWindow) || filterWindow < 0
%         filterWindow = 0;
%     end
% catch
%     filterWindow = 0;
% end
% 
% % Initialize StepChangesTable if empty / missing
% if ~isfield(app.Data, "StepChangesTable") || isempty(app.Data.StepChangesTable)
%     T = table( ...
%         double.empty(0,1), string.empty(0,1), string.empty(0,1), string.empty(0,1), ...
%         'VariableNames', {'Index','Axis','Action','String'} );
% else
%     T = app.Data.StepChangesTable;
%     % Normalize various possible stored types to a table with expected vars
%     if ~istable(T)
%         try
%             T = struct2table(T);
%         catch
%             try
%                 T = cell2table(T, 'VariableNames', {'Index','Axis','Action','String'});
%             catch
%                 T = table([],[],[],[], 'VariableNames', {'Index','Axis','Action','String'});
%             end
%         end
%     end
%     % Ensure expected variable names exist
%     missing = setdiff({'Index','Axis','Action','String'}, T.Properties.VariableNames);
%     for k = 1:numel(missing)
%         T.(missing{k}) = repmat("", height(T), 1);
%     end
% end
% 
% % Ensure Index numeric and Axis string
% if ~isnumeric(T.Index)
%     try
%         T.Index = double(T.Index);
%     catch
%         T.Index = nan(height(T),1);
%     end
% end
% T.Axis = string(T.Axis);
% 
% % Remove existing rows that are within filterWindow of any new stepIdx
% if ~isempty(T) && ~isempty(stepIdx)
%     existingIdx = T.Index(:);
%     existingAxis = T.Axis(:);
%     nExist = numel(existingIdx);
%     keepMask = true(nExist,1);
%     if filterWindow > 0
%         % For each existing row, check if there exists a new index within filterWindow
%         % and with the same axis label; if so, mark for removal.
%         for i = 1:nExist
%             if existingAxis(i) == axisVal
%                 if any(abs(existingIdx(i) - stepIdx) <= filterWindow)
%                     keepMask(i) = false;
%                 end
%             end
%         end
%     else
%         % filterWindow == 0: remove exact duplicates with same axis
%         sameAxisMask = (existingAxis == axisVal);
%         dupMask = ismember(existingIdx, stepIdx) & sameAxisMask;
%         keepMask(dupMask) = false;
%     end
%     T = T(keepMask, :);
% end
% 
% % Build and append new rows
% nNew = numel(stepIdx);
% Tnew = table(stepIdx, repmat(axisVal,nNew,1), repmat(actionVal,nNew,1), repmat(userStr,nNew,1), ...
%     'VariableNames', {'Index','Axis','Action','String'});
% 
% % Append and optionally sort
% T = [T; Tnew];
% T = sortrows(T, 'Index');
% 
% % Write back into app.Data.StepChangesTable
% app.Data.StepChangesTable = T;
% 
% setappdata(app.fig, 'app', app);
% 
% end
