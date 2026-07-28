function progFigure = parallelprogressdlg(varargin)
%parallelprogressdlg show/update progress dialog for a parfor/parallel job
%
% parallelprogressdlg creates or updates a progress dialog for a task running on
% parallel workers (e.g. inside a parfor loop). Because parallel workers cannot
% update the GUI directly, progress is inferred by periodically counting status
% files written by the workers (see 'Directory' below) using a timer, rather than
% by the caller calling parallelprogressdlg(n,...) on every iteration as with
% progressdlg.
%
% Like progressdlg, this is built entirely on documented, supported MATLAB
% graphics (uifigure/uiprogressdlg) rather than undocumented Java components.
% The timer periodically sets the ProgressDialog's Value/Message properties,
% which refresh automatically - no drawnow needed for those updates. A single
% drawnow is still issued right after the dialog is first created, to force its
% initial paint (see progressdlg.m for why).
%
% Syntax:
%   progressdlg                 creates default progress dialog
%
%   progressdlg('message');     creates default progress dialog with custom message
%
%   progressdlg('PropertyName',PropertyValue,...);  creates progress dialog with
%                                                   properties PropertyName/PropertyValue
%
%   progressdlg('close');       stops the refresh timer, removes the status
%                               directory, and closes the dialog
%
%   handle = progressdlg(...);  returns the handle to the progress dialog
%
%  PropertyName  | Description                                | Type           | Default
% ---------------|--------------------------------------------|----------------|------------------
%  Min           | Minimum value of ProgressBar               | double         | 0
%  Max           | Minimum value of ProgressBar               | double         | 100
%  Parent        | Parent Figure for Positioning              | figure_handle  | []
%  Position      | Override position of dialog                | [int x, int y] | []
%  ShowTime      | Displays an estimate of the remaining time | 'on/off'       | 'on'
%  Size          | Factor for dialog size (width=400px)       | double         | 1
%  String        | Descriptive text above ProgressBar         | char           | 'Progress'
%  Title         | Title of the dialog                        | char           | 'Progress'
%  Value         | Current Value of ProgressBar               | double         | 0
%  WindowStyle   | Normal or modal window behavior            | 'normal/modal' | 'modal'
%  Directory     | Filedirectory for saving progress          | char           | pwd
%
%
% Notes:
%   1) Requires MATLAB R2017b or later for uiprogressdlg.
%   2) The dialog is deliberately kept non-modal ('normal' WindowStyle) so the
%      main GUI remains usable while the parallel job runs - this matches the
%      previous version's behavior.
%   3) Performance is slightly decreased when using the remaining time estimation.
%
% Examples:
%   parallelprogressdlg('String','Analysing frames','Max',100,'Directory',pwd);
%   parfor n=1:100
%       % ... do the work, then signal progress by saving a status file ...
%       fSave(dirStatus,n);
%   end
%   parallelprogressdlg('close');
%
% See also: PROGRESSDLG, UIPROGRESSDLG, WAITBAR
%
% Bugs and suggestions:
%   Please send to Felix Ruhnow (felix.ruhnow at tu-dresden dot de)
%
% Adapted from:
% Daniel Claxton's WORKBAR
%
% Inspired by:
% Yair M. Altman's STATUSBAR
%
% Release history:
%    1.0  2013-06-10: initial version (undocumented Java components)
%    2.0  2026-07-20: rewritten on uifigure/uiprogressdlg for stability
%
% License to use and modify this code is granted freely to all interested, as long as the original author is
% referenced and attributed as such. The original author maintains the right to be solely associated with this work.
%
% Programmed and Copyright by Felix Ruhnow: felix.ruhnow(at)tu-dresden.de


    % keep variables for next call of progressdlg
    persistent p minv maxv uiHandles showtime;
    message = 'Progress';
    value = 0;
    text = [];

    % check input arguments
    if nargin < 1
        % progressdlg; creates default progress dialog
        uiHandles = [];

    elseif nargin == 1
        if ischar(varargin{1}) || (isstring(varargin{1}) && isscalar(varargin{1}))
            % progressdlg('close'); stops timer, cleans up status dir, closes dialog
            if strcmp(varargin{1},'close')
                closeDialog();
                progFigure = [];
                return
            else
                % progressdlg('message'); creates default progress dialog with custom message
                message = char(varargin{1});
                uiHandles = [];
            end

        elseif isnumeric(varargin{1})
            % progressdlg(n); creates default progress dialog with initial value n
            value = varargin{1};
        end
        varargin(:) = [];

    elseif nargin == 2 && isnumeric(varargin{1})
        % progressdlg(n,'message');
        if isempty(uiHandles)
            value = varargin{1};
            message = varargin{2};
            varargin(:) = [];
        else
            value = varargin{1};
            text = varargin{2}; %#ok<NASGU>
        end

    elseif isfield(uiHandles,'fig') && isvalid(uiHandles.fig)
        % progressdlg('PropertyName',PropertyValue,...); updates progress dialog
        parse(p,varargin{1:end});
        value = p.Results.value;
        text = p.Results.string; %#ok<NASGU>

    else
        % progressdlg('PropertyName',PropertyValue,...); creates progress dialog
        uiHandles = [];
    end

    if isempty(uiHandles)
        % if handles are empty create new dialog
        showtime = false;

        % get initial input arguments
        p = checkInput(value,message);
        parse(p,varargin{1:end});
        prog = p.Results;
        minv = prog.min;
        maxv = prog.max;

        if strcmp(prog.showtime,'on')
            showtime = true;
        end

        % size the host uifigure. uiprogressdlg has its own fixed internal layout
        % and is not scaled to fit a small host window, so give it a generous
        % minimum size (see progressdlg.m for the same reasoning)
        baseHeight = 200;
        if showtime
            baseHeight = baseHeight + 20;
        end
        width = ceil(max(400,400*prog.size));
        height = ceil(max(200,baseHeight*prog.size));

        % get position of progress dialog (middle of screen, if not otherwise)
        if isempty(prog.position)
            if isempty(prog.parent)
                screen = get(0,'ScreenSize');
                pos = ceil([(screen(3)-width)/2 (screen(4)-height)/2 width height]);
            else
                units = get(prog.parent,'Units');
                try
                    set(prog.parent,'Units','Pixels');
                catch
                    % uifigure parents always report Units in pixels and cannot be changed
                end
                pos_parent = get(prog.parent,'Position');
                try
                    set(prog.parent,'Units',units);
                catch
                end
                pos = ceil([pos_parent(1)+(pos_parent(3)-width)/2 pos_parent(2)+(pos_parent(4)-height)/2 width height]);
            end
        else
            pos = ceil([prog.position width height]);
        end

        % create the host uifigure that contains the progress dialog. WindowStyle
        % is deliberately left at its 'normal' default (not prog.windowstyle) so
        % the main GUI stays usable while the parallel job runs, matching the
        % previous version's behavior.
        fig = uifigure('Name',prog.title,'Position',pos,'Resize','off','Visible','off');
        fig.CloseRequestFcn = @(src,~) delete(src);

        % uiprogressdlg requires its parent uifigure to already be visible
        fig.Visible = 'on';

        msgLines = {prog.string};
        if showtime
            msgLines{2} = 'Estimated time remaining: ';
        end

        valueNorm = normalizeValue(prog.value,minv,maxv);

        dlg = uiprogressdlg(fig,'Title',prog.title,'Message',msgLines,...
                             'Indeterminate','off',...
                             'Cancelable','off',...
                             'ShowPercentage','on',...
                             'Value',valueNorm);

        % save handles of the created objects
        uiHandles.fig = fig;
        uiHandles.dlg = dlg;

        progFigure = fig;

        % force the dialog's first paint (see progressdlg.m for why this matters)
        drawnow;

        statusdir = [prog.directory filesep 'fiestastatus'];
        if isfolder(statusdir)
            rmdir(statusdir,'s');
        end
        mkdir(statusdir);

        % a plain (non-nested) callback is used deliberately: the timer fires
        % asynchronously, possibly long after this function has returned and its
        % persistent state may have been reused by another dialog, so all context
        % the callback needs is passed explicitly via UserData rather than
        % relying on a captured/shared workspace
        starttime = [];
        t = timer('TimerFcn',@RefreshProgressDlg,'Period',2,'TasksToExecute',Inf,'Tag','fParallel',...
                   'ExecutionMode','fixedRate','BusyMode','queue',...
                   'UserData',{starttime,minv,maxv,uiHandles,statusdir,showtime});
        start(t);
    else
        if value == maxv
            % if 'Maximum' value reached close and reset progress dialog
            closeDialog();
        end
    end

    function closeDialog()
        if ~isempty(uiHandles) && isfield(uiHandles,'dlg') && isvalid(uiHandles.dlg)
            close(uiHandles.dlg);
        end
        if ~isempty(uiHandles) && isfield(uiHandles,'fig') && isvalid(uiHandles.fig)
            delete(uiHandles.fig);
        end
        uiHandles = [];

        t = timerfind('Tag','fParallel');
        if ~isempty(t)
            var = get(t,'UserData');
            stop(t);
            delete(t);
            if iscell(var) && numel(var) >= 5 && isfolder(var{5})
                rmdir(var{5},'s');
            end
        end
    end

end

%% Clamp/normalize value into the [0,1] range expected by uiprogressdlg
function valueNorm = normalizeValue(value,minv,maxv)
    if maxv == minv
        valueNorm = 0;
    else
        valueNorm = (value-minv)/(maxv-minv);
    end
    valueNorm = min(1,max(0,valueNorm));
end

%% Timer callback: infer progress by counting status files written by workers
function RefreshProgressDlg(obj,~)
    var = get(obj,'UserData');
    starttime = var{1};
    minv = var{2};
    maxv = var{3};
    uiHandles = var{4};
    dirStatus = var{5};
    showtimeLocal = var{6};
    est_text = 'Estimated time remaining: ';

    if isempty(uiHandles) || ~isfield(uiHandles,'dlg') || ~isvalid(uiHandles.dlg)
        return
    end

    if isempty(starttime)
        starttime = clock;
        var{1} = starttime;
        set(obj,'UserData',var);
    end

    files = dir(fullfile(dirStatus,'*.mat'));
    value = length(files);

    if value < maxv
        try
            dlg = uiHandles.dlg;
            dlg.Value = normalizeValue(value,minv,maxv);

            % calculate and update remaining time if required
            if showtimeLocal && value > 0
                runtime = etime(clock,starttime);
                timeleft = runtime*(maxv-minv+1)/(value-minv+1) - runtime;
                timeleftstr = sec2timestr(timeleft);
                msg = dlg.Message;
                if ~iscell(msg)
                    msg = {msg};
                end
                msg{2} = [est_text timeleftstr];
                dlg.Message = msg;
            end
        catch
            return
        end
    end
end

%% Define input parameters and their default values
function p = checkInput(value,message)
    p = inputParser;
    states = {'on','off'};
    style = {'normal','modal'};
    addParameter(p,'title','Progress');
    addParameter(p,'size',1);
    addParameter(p,'value',value);
    addParameter(p,'min',0);
    addParameter(p,'max',100);
    addParameter(p,'position',[],@isnumeric);
    addParameter(p,'string',message);
    addParameter(p,'parent',[],@ishandle);
    addParameter(p,'showtime','on',@(x) any(validatestring(x,states)));
    addParameter(p,'windowstyle','modal',@(x) any(validatestring(x,style)));
    addParameter(p,'directory',pwd);
end

%% Convert seconds to hh:mm:ss
function timestr = sec2timestr(sec)
    h = floor(sec/3600); % Hours
    sec = sec - h*3600;
    m = floor(sec/60); % Minutes
    sec = sec - m*60;
    s = floor(sec); % Seconds
    timestr = sprintf('%02d:%02d:%02d',h,m,s);
end
