
function [x, y, button, ax] = ginputc(varargin)
%GINPUTC Graphical input from mouse.
%   GINPUTC behaves similarly to GINPUT, except you can customize the
%   cursor color, line width, and line style.
%
%   [X,Y] = GINPUTC(N) gets N points from the current axes and returns
%   the X- and Y-coordinates in length N vectors X and Y.  The cursor
%   can be positioned using a mouse.  Data points are entered by pressing
%   a mouse button or any key on the keyboard except carriage return,
%   which terminates the input before N points are entered.
%       Note: if there are multiple axes in the figure, use mouse clicks
%             instead of key presses. Key presses may not select the axes
%             where the cursor is.
%
%   [X,Y] = GINPUTC gathers an unlimited number of points until the return
%   key is pressed.
%
%   [X,Y] = GINPUTC(N, PARAM, VALUE) and [X,Y] = GINPUTC(PARAM, VALUE)
%   specifies additional parameters for customizing. Valid values for PARAM
%   are:
%       'FigHandle'     : Handle of the figure to activate. Default is gcf.
%       'Color'         : A three-element RGB vector, or one of the MATLAB
%                         predefined names, specifying the line color. See
%                         the ColorSpec reference page for more information
%                         on specifying color. Default is 'k' (black).
%       'LineWidth'     : A scalar number specifying the line width.
%                         Default is 0.5.
%       'LineStyle'     : '-', '--', '-.', ':'. Default is '-'.
%       'ShowPoints'    : TRUE or FALSE specifying whether to show the
%                         points being selected. Default is false.
%       'ConnectPoints' : TRUE or FALSE specifying whether to connect the
%                         points as they are being selected. This only
%                         applies when 'ShowPoints' is set to TRUE. Default
%                         is true.
%
%   [X,Y,BUTTON] = GINPUTC(...) returns a third result, BUTTON, that
%   contains a vector of integers specifying which mouse button was used
%   (1,2,3 from left) or ASCII numbers if a key on the keyboard was used.
%
%   [X,Y,BUTTON,AX] = GINPUTC(...) returns a fourth result, AX, that
%   contains a vector of axes handles for the data points collected.
%
%   Requires MATLAB R2007b or newer.
%
%   Examples:
%       [x, y] = ginputc;
%
%       [x, y] = ginputc(5, 'Color', 'r', 'LineWidth', 3);
%
%       [x, y, button] = ginputc(1, 'LineStyle', ':');
%
%       subplot(1, 2, 1); subplot(1, 2, 2);
%       [x, y, button, ax] = ginputc;
%
%       [x, y] = ginputc('ShowPoints', true, 'ConnectPoints', true);
%
%   See also GINPUT, GTEXT, WAITFORBUTTONPRESS.

% Jiro Doke
% October 19, 2012
% Copyright 2012 The MathWorks, Inc.

try
    if verLessThan('matlab', '7.5')
        error('ginputc:Init:IncompatibleMATLAB', ...
            'GINPUTC requires MATLAB R2007b or newer');
    end
catch %#ok<CTCH>
    error('ginputc:Init:IncompatibleMATLAB', ...
        'GINPUTC requires MATLAB R2007b or newer');
end

% Check input arguments
p = inputParser();

addOptional(p, 'N', inf, @(x) validateattributes(x, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}));
addParamValue(p, 'FigHandle', [], @(x) numel(x)==1 && ishandle(x));
addParamValue(p, 'Color', 'k', @colorValidFcn);
addParamValue(p, 'LineWidth', 0.5 , @(x) validateattributes(x, ...
    {'numeric'}, {'scalar', 'positive'}));
addParamValue(p, 'LineStyle', '-' , @(x) ischar(validatestring(x, ...
{'-', '--', '-.', ':'}))); 
addParamValue(p, 'ShowPoints', false, @(x) validateattributes(x, ...
    {'logical'}, {'scalar'}));
addParamValue(p, 'ConnectPoints', true, @(x) validateattributes(x, ...
    {'logical'}, {'scalar'}));

parse(p, varargin{:});

N = p.Results.N;
hFig = p.Results.FigHandle;
color = p.Results.Color;
linewidth = p.Results.LineWidth;
linestyle = p.Results.LineStyle;
showpoints = p.Results.ShowPoints;
connectpoints = p.Results.ConnectPoints;

%--------------------------------------------------------------------------
    function tf = colorValidFcn(in)
        % This function validates the color input parameter
        
        validateattributes(in, {'char', 'double'}, {'nonempty'});
        if ischar(in)
            validatestring(in, {'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'});
        else
            assert(isequal(size(in), [1 3]) && all(in>=0 & in<=1), ...
                'ginputc:InvalidColorValues', ...
                'RGB values for "Color" must be a 1x3 vector between 0 and 1');
            % validateattributes(in, {'numeric'}, {'size', [1 3], '>=', 0, '<=', 1})
        end
        tf = true;
    end
%--------------------------------------------------------------------------

if isempty(hFig)
    hFig = gcf;
end

% Try to get the current axes even if it has a hidden handle.
hAx = get(hFig, 'CurrentAxes');
if isempty(hAx)
    allAx = findall(hFig, 'Type', 'axes');
    if ~isempty(allAx)
        hAx = allAx(1);
    else
        hAx = axes('Parent', hFig);
    end
end

% Handle interactive properites of HG objects. Save the current settings so
% that they can be restored later
allHG = findall(hFig);
propsToChange = {...
    'WindowButtonUpFcn', ...
    'WindowButtonDownFcn', ...
    'WindowButtonMotionFcn', ...
    'WindowKeyPressFcn', ...
    'WindowKeyReleaseFcn', ...
    'ButtonDownFcn', ...
    'KeyPressFcn', ...
    'KeyReleaseFcn', ...
    'ResizeFcn'};
validObjects = false(length(allHG), length(propsToChange));
curCallbacks = cell(1, length(propsToChange));

% Save current properties and set them to ''
for id = 1:length(propsToChange)
    validObjects(:, id) = isprop(allHG, propsToChange{id});
    curCallbacks{id} = get(allHG(validObjects(:, id)), propsToChange(id));
    set(allHG(validObjects(:, id)), propsToChange{id}, '');
end

% Save current pointer
curPointer = get(hFig, 'Pointer');
curPointerShapeCData = get(hFig, 'PointerShapeCData');

% Change window functions
set(hFig, ...
    'WindowButtonDownFcn', @mouseClickFcn, ...
    'WindowButtonMotionFcn', @mouseMoveFcn, ...
    'KeyPressFcn', @keyPressFcn, ...
    'ResizeFcn', @resizeFcn, ...
    'Pointer', 'custom', ...
    'PointerShapeCData', nan(16, 16));

% Create an invisible axes for displaying the full crosshair cursor
hInvisibleAxes = axes(...
    'Parent', hFig, ...
    'Units', 'normalized', ...
    'Position', [0 0 1 1], ...
    'XLim', [0 1], ...
    'YLim', [0 1], ...
    'HitTest', 'off', ...
    'HandleVisibility', 'off', ...
    'Visible', 'off');

% Create line object for the selected points
if showpoints
    if connectpoints
        pointsLineStyle = '-';
    else
        pointsLineStyle = 'none';
    end
    
    selectedPoints = [];
    hPoints = line(nan, nan, ...
        'Parent', hInvisibleAxes, ...
        'HandleVisibility', 'off', ...
        'HitTest', 'off', ...
        'Color', [1 0 0], ...
        'Marker', 'o', ...
        'MarkerFaceColor', [1 .7 .7], ...
        'MarkerEdgeColor', [1 0 0], ...
        'LineStyle', pointsLineStyle);
end

% Create tooltip for displaying selected points
% hTooltipControl = text(0, 1, 'HIDE', ...
%     'Parent', hInvisibleAxes, ...
%     'HandleVisibility', 'callback', ...
%     'FontName', 'FixedWidth', ...
%     'VerticalAlignment', 'top', ...
%     'HorizontalAlignment', 'left', ...
%     'BackgroundColor', [.5 1 .5]);
% hTooltip = text(0, 0, 'No points', ...
%     'Parent', hInvisibleAxes, ...
%     'HandleVisibility', 'off', ...
%     'HitTest', 'off', ...
%     'FontName', 'FixedWidth', ...
%     'VerticalAlignment', 'top', ...
%     'HorizontalAlignment', 'left', ...
%     'BackgroundColor', [1 1 .5]);

% Call resizeFcn to update tooltip location
resizeFcn();

% Create full crosshair lines
hCursor = line(nan, nan, ...
    'Parent', hInvisibleAxes, ...
    'Color', color, ...
    'LineWidth', linewidth, ...
    'LineStyle', linestyle, ...
    'HandleVisibility', 'off', ...
    'HitTest', 'off');

% Prepare results
x = [];
y = [];
button = [];
ax = [];

% Wait until enter is pressed.
uiwait(hFig);


%--------------------------------------------------------------------------
    function mouseMoveFcn(varargin)
        % This function updates cursor location based on pointer location
        
        cursorPt = get(hInvisibleAxes, 'CurrentPoint');
        
        set(hCursor, ...
            'XData', [0 1 nan cursorPt(1) cursorPt(1)], ...
            'YData', [cursorPt(3) cursorPt(3) nan 0 1]);
    end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
    function mouseClickFcn(varargin)
        % This function captures mouse clicks.
        % If the tooltip control is clicked, then toggle tooltip display.
        % If anywhere else is clicked, record point.

%         if isequal(gco, hTooltipControl)
%             tooltipClickFcn();
%         else
            updatePoints(get(hFig, 'SelectionType'));
%         end
    end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
    function keyPressFcn(obj, edata) %#ok<INUSL>
        % This function captures key presses.
        % If "return", then exit.
        % If "delete" (or "backspace"), then delete previous point.
        % If any other key, record point.
        
        key = double(edata.Character);
        if isempty(key)
            return;
        end
        
        switch key
            case 13  % return
                exitFcn();
                
            case {8, 127}   % delete or backspace
                if ~isempty(x)
                    x(end) = [];
                    y(end) = [];
                    button(end) = [];
                    ax(end) = [];
                    
                    if showpoints
                        selectedPoints(end, :) = [];
                        set(hPoints, ...
                            'XData', selectedPoints(:, 1), ...
                            'YData', selectedPoints(:, 2));
                    end
                    
%                     displayCoordinates();
                end
                
            otherwise
                updatePoints(key);
                
        end
    end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
    function updatePoints(clickType)
        % This function captures the information for the selected point
        
        hAx = gca;
        pt = get(hAx, 'CurrentPoint');
        x = [x; pt(1)];
        y = [y; pt(3)];
        ax = [ax; hAx];

        if ischar(clickType)   % Mouse click
            switch lower(clickType)
                case 'open'
                    clickType = 1;
                case 'normal'
                    clickType = 1;
                case 'extend'
                    clickType = 2;
                case 'alt'
                    clickType = 3;
            end
        end
        button = [button; clickType];
        
%         displayCoordinates();
        
        if showpoints
            cursorPt = get(hInvisibleAxes, 'CurrentPoint');
            selectedPoints = [selectedPoints; cursorPt([1 3])];
            set(hPoints, ...
                'XData', selectedPoints(:, 1), ...
                'YData', selectedPoints(:, 2));
        end
        
        % If captured all points, exit
        if length(x) == N
            exitFcn();
        end
    end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%     function tooltipClickFcn()
%         % This function toggles the display of the tooltip
%         
%         if strcmp(get(hTooltipControl, 'String'), 'SHOW')
%             set(hTooltipControl, 'String', 'HIDE');
%             set(hTooltip, 'Visible', 'on');
%         else
%             set(hTooltipControl, 'String', 'SHOW');
%             set(hTooltip, 'Visible', 'off');
%         end
%     end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%     function displayCoordinates()
%         % This function updates the coordinates display in the tooltip
%         
%         if isempty(x)
%             str = 'No points';
%         else
%             str = sprintf('%d: %0.3f, %0.3f\n', [1:length(x); x'; y']);
%             str(end) = '';
%         end
%         set(hTooltip, ...
%             'String', str);
%     end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
    function resizeFcn(varargin)
        % This function adjusts the position of tooltip when the figure is
        % resized
        
%         sz = get(hTooltipControl, 'Extent');
%         set(hTooltip, 'Position', [0 sz(2)]);
    end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
    function exitFcn()
        % This function exits GINPUTC and restores previous figure settings
        
        for idx = 1:length(propsToChange)
            set(allHG(validObjects(:, idx)), propsToChange(idx), curCallbacks{idx});
        end
        
        % Restore window functions and pointer
        %         set(hFig, 'WindowButtonDownFcn', curWBDF);
        %         set(hFig, 'WindowButtonMotionFcn', curWBMF);
        %         set(hFig, 'WindowButtonUpFcn', curWBUF);
        %         set(hFig, 'KeyPressFcn', curKPF);
        %         set(hFig, 'KeyReleaseFcn', curKRF);
        %         set(hFig, 'ResizeFcn', curRF);
        
        % Restore pointer
        set(hFig, 'Pointer', curPointer);
        set(hFig, 'PointerShapeCData', curPointerShapeCData);

        % Delete invisible axes and return control
        delete(hInvisibleAxes);
        uiresume(hFig);
    end
%--------------------------------------------------------------------------

end