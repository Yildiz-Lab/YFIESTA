function options = fCheckboxGUI_MINFLUX(hMainGui, MFData)
% options = fCheckboxGUI_MINFLUX(hMainGui, MFData)
% Returns:
%   options.itr          - numeric vector of selected itr values
%   options.channelSplit - cell array of selected splits: {} or any subset of {'dcr','thi'}
%   options.dcrThreshold - numeric threshold for dcr (default 0.45)

% UI params
figW = 380; figH = 360; fontSize = 16; margin = 12; panelSpacing = 10; btnH = 44;
defaultDcr = 0.45;

% itr values (MFData optional)
if nargin < 2 || isempty(MFData)
    maxItr = 4;
else
    maxItr = max([MFData.itr]);
end
itrVals = 0:maxItr;
n = numel(itrVals);

% figure position (center on hMainGui if provided)
if nargin >= 1 && ishghandle(hMainGui)
    mpos = hMainGui.Position;
    x = round(mpos(1) + (mpos(3)-figW)/2);
    y = round(mpos(2) + (mpos(4)-figH)/2);
else
    scr = get(0,'ScreenSize');
    x = round((scr(3)-figW)/2);
    y = round((scr(4)-figH)/2);
end

fig = uifigure('Name','MINFLUX import parameters','Position',[x y figW figH]);
try theme(fig,"dark"); catch, end

% compute panel heights
availableH = figH - 2*margin - btnH - panelSpacing;
panelH1 = round(0.62 * availableH);
panelH2 = availableH - panelH1;

% Top panel: itr checkboxes
p1Y = figH - margin - panelH1;
p1 = uipanel(fig, 'Title','itr to import', 'Position',[margin p1Y figW-2*margin panelH1], 'FontSize', fontSize);

leftPad = 8; rightPad = 8; itemW = 64; itemH = 30; hSpacing = 8; vSpacing = 8;
x0 = leftPad;
y0 = p1.Position(4) - 40 - itemH;
maxRowWidth = p1.Position(3) - leftPad - rightPad;
cols = floor((maxRowWidth + hSpacing) / (itemW + hSpacing));
if cols < 1, cols = 1; end

chk = gobjects(n,1);
for k = 1:n
    col = mod(k-1, cols);
    row = floor((k-1)/cols);
    px = x0 + col*(itemW + hSpacing);
    py = y0 - row*(itemH + vSpacing);
    chk(k) = uicheckbox(p1, 'Text', num2str(itrVals(k)), ...
        'Position', [px py itemW itemH], 'FontSize', max(8,fontSize-4));
end
%set the last iteration to be accepted by default
chk(n).Value = true;

% Middle panel: channel splitting (multi-select checkboxes) + dcr threshold
p2Y = p1Y - panelSpacing - panelH2;
p2 = uipanel(fig, 'Title','channel splitting', 'Position',[margin p2Y figW-2*margin panelH2], 'FontSize', fontSize);

% Move these controls slightly lower by decreasing their y coordinates.
% Adjust the vertical offset here (increase moveDown to move further down).
moveDown = 12; % pixels to move controls down

% Multi-select checkboxes for channel splitting (default none selected)
cbDcr = uicheckbox(p2, 'Text', 'Split channels by dcr', ...
    'Position',[12 p2.Position(4)-44-moveDown 220 24], 'FontSize', max(8,fontSize-4), 'Value', false);
cbThi = uicheckbox(p2, 'Text', 'Split channels by thi', ...
    'Position',[12 p2.Position(4)-76-moveDown 220 24], 'FontSize', max(8,fontSize-4), 'Value', true);

% dcr threshold label + numeric edit field (placed to the right)
labelW = 110; editW = 80; ctrlH = 28;
lx = p2.Position(3) - rightPad - labelW - editW - 12;
ly = p2.Position(4) - 44 - moveDown;
uilabel(p2, 'Text', 'dcr threshold:', 'Position',[lx ly labelW ctrlH], 'FontSize', max(8,fontSize-6), 'HorizontalAlignment','right');
dcrEdit = uieditfield(p2, 'numeric', 'Value', defaultDcr, 'Limits', [0 inf], ...
    'Position', [lx+labelW+8 ly editW ctrlH], 'FontSize', max(8,fontSize-6));

% Apply button
uibutton(fig, 'Text','Apply', 'FontSize', fontSize, ...
    'Position',[margin 8 figW-2*margin btnH], 'ButtonPushedFcn',@applyCb);

% default output if closed without Apply
options = struct('itr',[], 'channelSplit',{{}}, 'dcrThreshold', defaultDcr);

% block until Apply or figure closed
uiwait(fig);

% after resume: read stored options if set
if isvalid(fig) && isprop(fig,'UserData') && isstruct(fig.UserData)
    options = fig.UserData;
end
if isvalid(fig)
    close(fig);
end

    function applyCb(~,~)
        % collect selected itr values
        selected = itrVals([chk.Value]);
        % collect channel split selections into cell array (possibly empty)
        sel = {};
        if cbDcr.Value, sel{end+1} = 'dcr'; end
        if cbThi.Value, sel{end+1} = 'thi'; end
        % read dcr threshold (validate numeric)
        thr = defaultDcr;
        try
            val = dcrEdit.Value;
            if isnumeric(val) && isscalar(val) && ~isnan(val)
                thr = double(val);
            end
        catch
            thr = defaultDcr;
        end
        fig.UserData = struct('itr', selected, 'channelSplit', {sel}, 'dcrThreshold', thr);
        uiresume(fig);
    end
end

