function app = UserStepChangesTable(app, stepIdx, mode, varargin)
% UserStepChangesTable  Maintain a state table of currently-saved step changes.
%
% The table reflects which steps ARE saved (not which actions were taken).
% Deletions remove rows; additions insert rows. The Action column is gone.
%
% Table schema: Index (double) | Axis (string) | String (string)
%
% Modes:
%   "fit"   - Add detected steps (String = "fit"). Existing rows for the
%             same Axis within filterWindow are replaced.
%   "user"  - Add or delete a user-specified step depending on BtnAdd /
%             BtnDelete state (or explicit "Action" name-value arg).
%             "add"    -> insert row with provided UserString
%             "delete" -> remove matching row(s)
%   "align" - Shift existing step indices after an alignment operation.
%             Rows within filterWindow of a given stepIdx get their Index
%             updated to that new value; String is preserved.
%             If actionVal == "delete", matched rows are removed instead.
%
% Name-Value inputs:
%   "Action"     - "add" | "delete"  (overrides button inference)
%   "UserString" - string stored in the String column (default: "")

% ── Guard ────────────────────────────────────────────────────────────────
if isempty(stepIdx)
    return
end
stepIdx = double(stepIdx(:));   % ensure column vector of doubles

% ── Parse name-value args ─────────────────────────────────────────────────
p = inputParser;
addParameter(p, "Action",     "", @(x) ischar(x) || isstring(x));
addParameter(p, "UserString", "", @(x) ischar(x) || isstring(x));
parse(p, varargin{:});
actionVal = string(p.Results.Action);
userStr   = string(p.Results.UserString);

% ── Axis ──────────────────────────────────────────────────────────────────
axisVal = "unknown";
try
    axisVal = string(app.FitPanelFields.DropFitAxis.Value);
catch
end

% ── Resolve mode-specific defaults ───────────────────────────────────────
switch string(mode)
    case "fit"
        actionVal = "add";
        userStr   = "fit";

    case {"user", "align"}
        % Infer action from buttons if not supplied (used by "user" mode only)
        if actionVal == ""
            try
                if app.FitPanelFields.BtnDelete.UserData.pressed
                    actionVal = "delete";
                end
            catch
            end
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
        if userStr == "" && actionVal ~= "delete"
            try
                userStr = string(app.FitPanelFields.DropFitString.Value);
            catch
            end
        end

    otherwise
        error("Mode must be ''fit'', ''user'', or ''align''.");
end

% ── Filter window ─────────────────────────────────────────────────────────
filterWindow = 0;
try
    fw = double(app.FilterPanelFields.FilterWindow.Value);
    if isfinite(fw) && fw >= 0
        filterWindow = fw;
    end
catch
end

% ── Load / initialize state table ────────────────────────────────────────
emptyT = table( ...
    double.empty(0,1), string.empty(0,1), string.empty(0,1), ...
    'VariableNames', {'Index','Axis','String'} );

if ~isfield(app.Data, "StepChangesTable") || isempty(app.Data.StepChangesTable)
    T = emptyT;
else
    T = app.Data.StepChangesTable;
    % Coerce to table if stored as struct or cell
    if ~istable(T)
        try,   T = struct2table(T);
        catch
            try, T = cell2table(T, 'VariableNames', {'Index','Axis','String'});
            catch
                T = emptyT;
            end
        end
    end
    % Drop legacy Action column if present
    if ismember('Action', T.Properties.VariableNames)
        T.Action = [];
    end
    % Ensure all expected columns exist
    for col = {'Index','Axis','String'}
        if ~ismember(col{1}, T.Properties.VariableNames)
            if strcmp(col{1}, 'Index')
                T.(col{1}) = nan(height(T), 1);
            else
                T.(col{1}) = repmat("", height(T), 1);
            end
        end
    end
end

% ── Normalize column types ────────────────────────────────────────────────
if ~isnumeric(T.Index)
    try,   T.Index = double(T.Index);
    catch, T.Index = nan(height(T), 1);
    end
end
T.Axis   = string(T.Axis);
T.String = string(T.String);

% ── Main logic per mode ───────────────────────────────────────────────────
switch string(mode)

    % ── FIT ──────────────────────────────────────────────────────────────
    case "fit"
        % Remove any existing rows on this Axis within filterWindow of new indices
        if ~isempty(T)
            keepMask = true(height(T), 1);
            for i = 1:height(T)
                if T.Axis(i) == axisVal
                    if filterWindow > 0
                        if any(abs(T.Index(i) - stepIdx) <= filterWindow)
                            keepMask(i) = false;
                        end
                    else
                        if ismember(T.Index(i), stepIdx)
                            keepMask(i) = false;
                        end
                    end
                end
            end
            T = T(keepMask, :);
        end
        % Append new rows
        nNew = numel(stepIdx);
        T = [T; table(stepIdx, repmat(axisVal,nNew,1), repmat(userStr,nNew,1), ...
            'VariableNames', {'Index','Axis','String'})];

    % ── USER ─────────────────────────────────────────────────────────────
    case "user"
        if actionVal == "delete"
            % Remove rows matching exact index on this Axis.
            % filterWindow is intentionally ignored here — the user is
            % targeting a specific step, not a neighborhood of steps.
            keepMask = ~(T.Axis == axisVal & ismember(T.Index, stepIdx));
            T = T(keepMask, :);
        else
            % "add": remove any existing row at the exact same index on this
            % Axis, then insert. filterWindow is not used here — the user is
            % placing a step at a specific index and should not evict neighbors.
            keepMask = ~(T.Axis == axisVal & ismember(T.Index, stepIdx));
            T = T(keepMask, :);
            nNew = numel(stepIdx);
            T = [T; table(stepIdx, repmat(axisVal,nNew,1), repmat(userStr,nNew,1), ...
                'VariableNames', {'Index','Axis','String'})];
        end

    % ── ALIGN ─────────────────────────────────────────────────────────────
    case "align"
        % Shift existing step indices to new positions only.
        % No rows are inserted or deleted — only Index values are updated.
        % String is preserved. Axis is updated to current axisVal.
        for k = 1:numel(stepIdx)
            newIdx = stepIdx(k);
            if filterWindow > 0
                matchMask = (T.Axis == axisVal) & (abs(T.Index - newIdx) <= filterWindow);
            else
                matchMask = (T.Axis == axisVal) & (T.Index == newIdx);
            end
            matches = find(matchMask);
            for m = 1:numel(matches)
                T.Index(matches(m)) = newIdx;
                T.Axis(matches(m))  = axisVal;
            end
        end

        % Collapse exact-duplicate rows that may result from multiple steps
        % shifting onto the same index (keep first occurrence)
        if height(T) > 1
            key = strcat(string(T.Index), "|", T.Axis, "|", T.String);
            [~, ia] = unique(key, 'stable');
            keepMask = false(height(T), 1);
            keepMask(ia) = true;
            T = T(keepMask, :);
        end
end

% ── Sort and write back ───────────────────────────────────────────────────
try
    T = sortrows(T, 'Index');
catch
end

% Defensive strip of Action column in case anything upstream reintroduced it
if ismember('Action', T.Properties.VariableNames)
    T.Action = [];
end

app.Data.StepChangesTable = T;
try
    app.StepChangesTable = T;
catch
end

% T

try
    setappdata(app.fig, 'app', app);
catch
end

end