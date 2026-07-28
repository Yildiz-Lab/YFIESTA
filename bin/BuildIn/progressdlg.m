function progFigure = progressdlg(varargin)
%progressdlg show/update progress dialog
%
% progressdlg creates or updates an existing progress dialog that reports on the
% progress of a task including the estimated time remaining as well as a button to
% cancel the task. It is built entirely on documented, supported MATLAB graphics
% (uifigure/uiprogressdlg), so unlike the previous version it does not need drawnow
% to refresh the display: uifigure-based graphics update automatically whenever a
% property is set.
%
% Syntax:
%   progressdlg                 creates default progress dialog
%
%   progressdlg('message');     creates default progress dialog with custom message
%
%   progressdlg(n);             creates default progress dialog with initial value n
%                               or updates existing progress dialog to value n
%
%   progressdlg(n,'message');   creates default progress dialog with message and value n
%                               or updates existing progress dialog to value n and message
%
%   progressdlg('PropertyName',PropertyValue,...);  creates or updates progress dialog with
%                                                   properties PropertyName/PropertyValue
%
%   handle = progressdlg(...);  returns the handle to the progress dialog
%                               handle is empty if user aborted by pressing cancel
%
%  PropertyName  | Description                                | Type           | Default
% ---------------|--------------------------------------------|----------------|------------------
%  Cancel        | Displays a cancel button to abort          | 'on,off'       | 'off'
%  Indeterminate | Sets the ProgressBar to indeterminate      | 'on,off'       | 'off'
%  Min           | Minimum value of ProgressBar               | double         | 0
%  Max           | Minimum value of ProgressBar               | double         | 100
%  Parent        | Parent Figure for Positioning              | figure_handle  | []
%  Position      | Override position of dialog                | [int x, int y] | []
%  ShowTime      | Displays an estimate of the remaining time | 'on/off'       | 'on'
%  Size          | Factor for dialog size (width=400px)       | double         | 1
%  String        | Descriptive text above ProgressBar         | char           | 'Progress'
%  Title         | Title of the dialog                        | char           | 'Progress'
%  Value         | Current Value of ProgressBar               | double         | 0
%  WindowStyle   | Normal or modal window behavior             | 'normal/modal' | 'modal'
%
%
% Notes:
%   1) Requires MATLAB R2017b or later for uiprogressdlg. The 'modal' WindowStyle
%      option requires R2020b or later; on earlier releases it silently falls back
%      to 'normal'.
%   2) The host uifigure is created solely to contain the built-in progress dialog,
%      mirroring the way the previous version created a dedicated figure.
%
% Examples:
%   progressdlg('Calculating magic square');
%   for n=1:100
%       magic(1000);
%       progressdlg(n);
%   end
%
%   progressdlg('Title','Progress Dialog','String','Calculating magic squares','Min',0,'Max',1000,'Value',0);
%   for n=1:1000
%       magic(400);
%       progressdlg(n);
%   end
%
%   progressdlg('Title','Progress Dialog','String','Calculating magic squares','Position',[100 100],'WindowStyle','normal');
%   for n=1:10
%       magic(10000);
%       progressdlg(n*10,['Calculating magic square - Round: ' num2str(n)]);
%   end
%
%   h = progressdlg('Title','Progress','String','Calculating magic square','ShowTime','on','Cancel','on','Size',1.5);
%   for n=1:100
%       magic(4000);
%       h = progressdlg(n);
%       if isempty(h)
%           break;
%       end
%   end
%
% See also: UIPROGRESSDLG, WAITBAR
%
% Adapted from:
% Daniel Claxton's WORKBAR
%
% Inspired by:
% Yair M. Altman's STATUSBAR
%
% Release history:
%    1.0  2013-06-10: initial version (undocumented Java components)
%    2.0  2026-07-17: rewritten on uifigure/uiprogressdlg for stability
%
% License to use and modify this code is granted freely to all interested, as long as the original author is
% referenced and attributed as such. The original author maintains the right to be solely associated with this work.
%
% Programmed and Copyright by Felix Ruhnow: felix.ruhnow(at)tu-dresden.de


    % keep variables for next call of progressdlg
    persistent p minv maxv cancelOn uiHandles showtime starttime msgLines;
    message = 'Progress';
    value = 0;
    text = [];
    est_text = 'Estimated time remaining: ';

    % check input arguments
    if nargin < 1
        % progressdlg; creates default progress dialog
        uiHandles = [];

    elseif nargin == 1
        if ischar(varargin{1}) || (isstring(varargin{1}) && isscalar(varargin{1}))
            % progressdlg('close'); closes progress dialog
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
            %                 or updates existing progress dialog to value n
            value = varargin{1};
        end
        varargin(:) = [];

    elseif nargin == 2 && isnumeric(varargin{1})
        % progressdlg(n,'message');
        if isempty(uiHandles)
            % creates default progress dialog with message and value n
            value = varargin{1};
            message = varargin{2};
            varargin(:) = [];
        else
            % updates existing progress dialog to value n and message
            value = varargin{1};
            text = varargin{2};
        end

    elseif isfield(uiHandles,'fig') && isvalid(uiHandles.fig)
        % progressdlg('PropertyName',PropertyValue,...); updates progress dialog
        parse(p,varargin{1:end});
        value = p.Results.value;
        text = p.Results.string;

    else
        % progressdlg('PropertyName',PropertyValue,...); creates progress dialog
        uiHandles = [];
    end

    if isempty(uiHandles)
        % if handles are empty create new dialog
        cancelOn = false;
        showtime = false;

        % get initial input arguments
        p = checkInput(value,message);
        parse(p,varargin{1:end});
        prog = p.Results;
        minv = prog.min;
        maxv = prog.max;

        if strcmp(prog.showtime,'on') && ~strcmp(prog.indeterminate,'on')
            showtime = true;
        end
        if strcmp(prog.cancel,'on')
            cancelOn = true;
        end

        % size the host uifigure. uiprogressdlg has its own fixed internal layout
        % (title, message lines, the progress bar itself, and the cancel button),
        % and it is NOT scaled to fit a small host window - if the host figure is
        % too small, the bar gets clipped and doesn't render. So we give it a
        % generous minimum size and only grow it for extra message lines / cancel.
        baseHeight = 200;
        if showtime
            baseHeight = baseHeight + 20;
        end
        if cancelOn
            baseHeight = baseHeight + 30;
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

        % create the host uifigure that contains the progress dialog
        fig = uifigure('Name',prog.title,'Position',pos,'Resize','off','Visible','off');
        try
            fig.WindowStyle = prog.windowstyle;
        catch
            % 'modal' WindowStyle for uifigure requires R2020b or later; ignore on older releases
        end
        fig.CloseRequestFcn = @(src,~) delete(src);

        % uiprogressdlg requires its parent uifigure to already be visible
        fig.Visible = 'on';

        % build the (optionally two-line) message: descriptive text + estimated time remaining
        msgLines = {prog.string};
        if showtime
            msgLines{2} = est_text;
        end

        showPercentage = 'on';
        if strcmp(prog.indeterminate,'on')
            showPercentage = 'off';
        end

        valueNorm = normalizeValue(prog.value,minv,maxv);

        dlg = uiprogressdlg(fig,'Title',prog.title,'Message',msgLines,...
                             'Indeterminate',prog.indeterminate,...
                             'Cancelable',prog.cancel,...
                             'ShowPercentage',showPercentage,...
                             'Value',valueNorm);

        % save handles of the created objects
        uiHandles.fig = fig;
        uiHandles.dlg = dlg;

        progFigure = fig;

        % force the dialog's first paint. Without this, if the calling code goes
        % straight into a tight loop with no other graphics/event activity, the
        % dialog can remain visually blank (this mirrors MathWorks' own
        % uiprogressdlg examples, which always call drawnow right after creation)
        drawnow;

        % save start time for progress
        starttime = clock;
    else
        % if handles are present update progress dialog
        fig = uiHandles.fig;
        dlg = uiHandles.dlg;

        if ~isvalid(fig) || ~isvalid(dlg)
            % dialog was closed some other way (e.g. window close button)
            uiHandles = [];
            progFigure = [];
            return
        end

        % check for cancel button before doing any further work
        if cancelOn && dlg.CancelRequested
            closeDialog();
            progFigure = [];
            return
        end

        % update ProgressBar
        dlg.Value = normalizeValue(value,minv,maxv);

        % update descriptive text if required
        if ~isempty(text)
            msgLines{1} = text;
        end

        % calculate and update remaining time if required
        if showtime
            if value == minv
                starttime = clock;
            else
                runtime = etime(clock,starttime);
                timeleft = runtime*(maxv-minv)/(value-minv) - runtime;
                timeleftstr = sec2timestr(timeleft);
                msgLines{2} = [est_text timeleftstr];
            end
        end

        if ~isempty(text) || showtime
            dlg.Message = msgLines;
        end

        % check again for cancel button after the update
        if cancelOn && dlg.CancelRequested
            closeDialog();
            progFigure = [];
            return
        end

        progFigure = fig;
    end

    if value == maxv
        % if 'Maximum' value reached close and reset progress dialog
        closeDialog();
        progFigure = [];
    end

    function closeDialog()
        if ~isempty(uiHandles) && isfield(uiHandles,'dlg') && isvalid(uiHandles.dlg)
            close(uiHandles.dlg);
        end
        if ~isempty(uiHandles) && isfield(uiHandles,'fig') && isvalid(uiHandles.fig)
            delete(uiHandles.fig);
        end
        uiHandles = [];
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

%% Define input parameters and their default values
function p = checkInput(value,message)
    p = inputParser;
    states = {'on','off'};
    style = {'normal','modal'};
    addParamValue(p,'title','Progress');
    addParamValue(p,'size',1);
    addParamValue(p,'value',value);
    addParamValue(p,'min',0);
    addParamValue(p,'max',100);
    addParamValue(p,'position',[],@isnumeric);
    addParamValue(p,'string',message);
    addParamValue(p,'parent',[],@ishandle);
    addParamValue(p,'cancel','off',@(x) any(validatestring(x,states)))
    addParamValue(p,'showtime','on',@(x) any(validatestring(x,states)))
    addParamValue(p,'indeterminate','off',@(x) any(validatestring(x,states)))
    addParamValue(p,'windowstyle','normal',@(x) any(validatestring(x,style)))
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