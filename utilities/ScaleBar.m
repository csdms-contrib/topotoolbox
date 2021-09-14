classdef ScaleBar < handle
    
%ScaleBar A dynamic scalebar
%
% Syntax
%
%     SB = ScaleBar
%     SB = ScaleBar(pn,pv,...)
%
% Description
%
%     ScaleBar plots a dynamic scale bar in the 2D-axis. The function is
%     currently beta and only supports axis with projected coordinate
%     systems (no geographic coordinates). ScaleBar refers to the x-axis of
%     a plot and may be used for plots with different axis aspect ratios.
%
%     Once the ScaleBar SB is created, you can make changes by directly
%     setting the properties, e.g.
%
%     SB = ScaleBar;
%     SB.location = 'northwest';
%
%     In addition, there are functions to increase or decrease the font
%     size.
%
%     smaller(SB)
%     bigger(SB)
%     
% Input arguments
%
%     'ax'           axes handle {gca}
%     'color'        color of text and scalebar {'k'}
%     'xyunit'       unit of coordinates {'m'} or 'km' or 'none'
%     'displayunit'  unit shown by scale bar label {'auto'}, 'm', 'km' or
%                    'none'. If set to 'auto', then the label switch
%                    automatically between 'm' and 'km'.
%     'rellength'    length of the scale bar relative to the length of the
%                    x axis. This property controls the length of the scale
%                    bar. It should range between 0.1 and 0.5.
%     'location'     placement of scalebar inside the axes {'southeast'},
%                    'southwest','northeast', or 'northwest'.
%     'backgroundcolor'   background color of the text box {'none'}.
%     'fontsize'     Fontsize {10}
%
% Output arguments
%
%     SB     handle to scale bar
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     imageschs(DEM)
%     SB = ScaleBar;
%     SB.location = 'northeast';
%     SB.color = 'w';
%     smaller(SB,2)
%
% References: The function uses plotboxpos by Kelly Kearney. The function
% is available here: https://github.com/kakearney/plotboxpos-pkg
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 8. September, 2021

    properties(SetAccess = 'private')
        scale
        scaletext
        ax
        el
    end
    
    properties(GetAccess = 'public', SetAccess = 'public')
        color       = 'k'
        displayunit = 'auto'
        xyunit      = 'm'
        rellength   = 0.2;
        location    = 'southeast'
        backgroundcolor = 'none'
        fontsize    = 10;        
    end
    
    methods
        
        function SB = ScaleBar(varargin)
            
            % Parse inputs
            p = inputParser;
            addParameter(p,'parent',gca);
            addParameter(p,'xyunit','m');
            addParameter(p,'displayunit','auto');
            addParameter(p,'color','k');
            addParameter(p,'rellength',0.2);
            addParameter(p,'backgroundcolor','none');
            addParameter(p,'fontsize',10);
            addParameter(p,'location','southeast',@(x) ischar(validatestring(x,{'northwest','southwest','northeast','southeast'})));
            parse(p,varargin{:});
            
            params = p.Results;
            
            % set initial location of scalebar
            x1 = 0;
            x2 = 1;
            y1 = 0;
            y2 = 0;
            
            SB.scaletext = annotation('textbox','EdgeColor','none',...
                'BackgroundColor',params.backgroundcolor,...
                'Position',[y1 y2 1 0.05],...
                'String','1',...
                'HorizontalAlignment','center',...
                'VerticalAlignment','bottom',...
                'FitBoxToText','on',...
                'Fontsize',params.fontsize,...
                'Margin',3,...
                'Color',params.color);
            SB.scale     = annotation('line',[x1 x2], [y1 y2],'Color',params.color,'LineWidth',2);
            SB.ax     = params.parent;
            SB.xyunit = params.xyunit;
            SB.displayunit = params.displayunit;
            SB.location  = params.location;
            SB.rellength = params.rellength;
            
            %% Create listeners
            hfig = SB.ax.Parent;
            SB.el = addlistener(SB.ax,{'XLim','YLim', 'Position', 'OuterPosition'},'PostSet',@(src,evnt)placescalebar(SB));
            hfig.SizeChangedFcn = @(src,evnt) placescalebar(SB);    
        end
        
        function set.color(SB,c)
            colorize(SB,'color',c)
        end
        function set.backgroundcolor(SB,c)
            colorize(SB,'backgroundcolor',c)
        end
        function set.location(SB,location)
            location = validatestring(location,{'northwest','southwest','northeast','southeast'});
            SB.location = location;
            placescalebar(SB)
        end
        function set.rellength(SB,l)
            validateattributes(l,'numeric',{'>',0,'<',0.8})
            SB.rellength = l;
            placescalebar(SB)
        end
        function set.xyunit(SB,unit)
            unit = validatestring(unit,{'m','km','none'});
            SB.xyunit = unit;
            placescalebar(SB)
        end
        function set.displayunit(SB,unit)
            unit = validatestring(unit,{'m','km','none','auto'});
            SB.displayunit = unit;
            placescalebar(SB)
        end
        function set.fontsize(SB,fs)
            textchange(SB,'FontSize',fs);
        end
        
        function delete(SB)
            % Delete ScaleBar
            
            delete(SB.el)
            delete(SB.scale)
            delete(SB.scaletext)
            try
            hfig = SB.ax.Parent; 
            hfig.SizeChangedFcn = [];
            end
            
        end
        
        function bigger(SB,val)
            % Increase font size of scale bar text
            if nargin == 1
                textchange(SB,'Bigger',1)
            else
                textchange(SB,'Bigger',val)
            end
        end
        function smaller(SB,val)
            % Decrease font size of scale bar text
            if nargin == 1
                textchange(SB,'Smaller',1)
            else
                textchange(SB,'Smaller',val)
            end
        end
        
        function textchange(SB,type,val)
            switch type
                case 'FontSize'
                    SB.scaletext.FontSize = val;
                case 'Bigger'
                    SB.scaletext.FontSize = SB.scaletext.FontSize + val;
                case 'Smaller'
                    SB.scaletext.FontSize = SB.scaletext.FontSize - val;
            end
            placescalebar(SB)
        end
             
        function colorize(SB,type,c)
            switch type
                case 'color'
                    SB.scale.Color = c;
                    SB.scaletext.Color = c;
                case 'backgroundcolor'
                    SB.scaletext.BackgroundColor = c;
            end
        end
        
        
        function placescalebar(SB)
            
            lims = axis(SB.ax);
            
            units_former = SB.ax.Units;
            SB.ax.Units = 'normalized';
            pos  = plotboxpos(SB.ax);
            
            xext = lims(2)-lims(1);
            scalefactor = pos(3)/xext;
            sblengthnorm = SB.rellength; % scalebar should have about 1/10 of the x-axis extent
            sblength = xext*sblengthnorm;
            sblengthvalues = [0.1 0.25 .5 1] * 10.^(ceil(log10(sblength)));
            % find nearest value
            [~,ix] = min(abs(sblengthvalues-sblength));
            sblength = sblengthvalues(ix);
            
            switch SB.xyunit
                case 'none'
                    sblengthtext = num2str(sblength);
                case 'm'                   
                    switch SB.displayunit
                        case 'm'
                            sblengthtext = [num2str(sblength) ' m'];
                        case 'km'
                            sblengthtext = [num2str(sblength/1000) ' km'];
                        case 'auto'
                            if sblength > 1000
                                sblengthtext = [num2str(sblength/1000) ' km'];
                            else
                                sblengthtext = [num2str(sblength) ' m'];
                            end
                        case 'none'
                            sblengthtext = num2str(sblength);
                    end
                case 'km'
                    switch SB.displayunit
                        case 'm'
                            sblengthtext = [num2str(sblength*1000) ' m'];
                        case 'km'
                            sblengthtext = [num2str(sblength) ' km'];
                        case 'auto'
                            if sblength > 1000
                                sblengthtext = [num2str(sblength) ' km'];
                            else
                                sblengthtext = [num2str(sblength*1000) ' m'];
                            end
                        case 'none'
                            sblengthtext = num2str(sblength);
                    end
            end
            
            sblength = sblength*scalefactor;
            
            % nearest values
            switch SB.location
                case 'southeast'
                    x1 = pos(1) + pos(3)*0.95;
                    y1 = pos(2) + pos(4)*0.05;
                    x2 = x1-sblength;
                    y2 = y1;
                    x  = [x2 x1];
                    y  = [y1 y2];
                case 'northeast'
                    x1 = pos(1) + pos(3)*0.95;
                    y1 = pos(2) + pos(4)*0.85;
                    x2 = x1-sblength;
                    y2 = y1;
                    x  = [x2 x1];
                    y  = [y1 y2];
                case 'northwest'
                    x1 = pos(1) + pos(3)*0.05;
                    y1 = pos(2) + pos(4)*0.85;
                    x2 = x1+sblength;
                    y2 = y1;
                    x  = [x1 x2];
                    y  = [y1 y2];
                case 'southwest'
                    x1 = pos(1) + pos(3)*0.05;
                    y1 = pos(2) + pos(4)*0.05;
                    x2 = x1+sblength;
                    y2 = y1;
                    x  = [x1 x2];
                    y  = [y1 y2];
                    
            end
            
            SB.scale.X = [x1 x2];
            SB.scale.Y = [y1 y2];
            SB.scaletext.FitBoxToText = 'on';
            textheight = SB.scaletext.Position(4);
            SB.scaletext.Position = [min(x) y(1) sblength textheight];
            pr = SB.scaletext.Units;
            SB.scaletext.Units = 'Points';
            SB.scaletext.Position(4) = SB.scaletext.FontSize + SB.scaletext.Margin*2;
            SB.scaletext.Units = pr;
            SB.scaletext.String = sblengthtext;
            
            SB.ax.Units = units_former;
                       
        end
        
        
    end

end

function pos = plotboxpos(h)
%PLOTBOXPOS Returns the position of the plotted axis region
%
% pos = plotboxpos(h)
%
% This function returns the position of the plotted region of an axis,
% which may differ from the actual axis position, depending on the axis
% limits, data aspect ratio, and plot box aspect ratio.  The position is
% returned in the same units as the those used to define the axis itself.
% This function can only be used for a 2D plot.  
%
% Input variables:
%
%   h:      axis handle of a 2D axis (if ommitted, current axis is used).
%
% Output variables:
%
%   pos:    four-element position vector, in same units as h

% Copyright 2010 Kelly Kearney

% Check input

if nargin < 1
    h = gca;
end

if ~ishandle(h) || ~strcmp(get(h,'type'), 'axes')
    error('Input must be an axis handle');
end

% Get position of axis in pixels

currunit = get(h, 'units');
axisPos  = getpixelposition(h);

% Calculate box position based axis limits and aspect ratios

darismanual  = strcmpi(get(h, 'DataAspectRatioMode'),    'manual');
pbarismanual = strcmpi(get(h, 'PlotBoxAspectRatioMode'), 'manual');

if ~darismanual && ~pbarismanual
    
    pos = axisPos;
    
else

    xlim = get(h, 'XLim');
    ylim = get(h, 'YLim');
    
    % Deal with axis limits auto-set via Inf/-Inf use
    
    if any(isinf([xlim ylim]))
        hc = get(h, 'Children');
        hc(~arrayfun( @(h) isprop(h, 'XData' ) & isprop(h, 'YData' ), hc)) = [];
        xdata = get(hc, 'XData');
        if iscell(xdata)
            xdata = cellfun(@(x) x(:), xdata, 'uni', 0);
            xdata = cat(1, xdata{:});
        end
        ydata = get(hc, 'YData');
        if iscell(ydata)
            ydata = cellfun(@(x) x(:), ydata, 'uni', 0);
            ydata = cat(1, ydata{:});
        end
        isplotted = ~isinf(xdata) & ~isnan(xdata) & ...
                    ~isinf(ydata) & ~isnan(ydata);
        xdata = xdata(isplotted);
        ydata = ydata(isplotted);
        if isempty(xdata)
            xdata = [0 1];
        end
        if isempty(ydata)
            ydata = [0 1];
        end
        if isinf(xlim(1))
            xlim(1) = min(xdata);
        end
        if isinf(xlim(2))
            xlim(2) = max(xdata);
        end
        if isinf(ylim(1))
            ylim(1) = min(ydata);
        end
        if isinf(ylim(2))
            ylim(2) = max(ydata);
        end
    end

    dx = diff(xlim);
    dy = diff(ylim);
    dar = get(h, 'DataAspectRatio');
    pbar = get(h, 'PlotBoxAspectRatio');

    limDarRatio = (dx/dar(1))/(dy/dar(2));
    pbarRatio = pbar(1)/pbar(2);
    axisRatio = axisPos(3)/axisPos(4);

    if darismanual
        if limDarRatio > axisRatio
            pos(1) = axisPos(1);
            pos(3) = axisPos(3);
            pos(4) = axisPos(3)/limDarRatio;
            pos(2) = (axisPos(4) - pos(4))/2 + axisPos(2);
        else
            pos(2) = axisPos(2);
            pos(4) = axisPos(4);
            pos(3) = axisPos(4) * limDarRatio;
            pos(1) = (axisPos(3) - pos(3))/2 + axisPos(1);
        end
    elseif pbarismanual
        if pbarRatio > axisRatio
            pos(1) = axisPos(1);
            pos(3) = axisPos(3);
            pos(4) = axisPos(3)/pbarRatio;
            pos(2) = (axisPos(4) - pos(4))/2 + axisPos(2);
        else
            pos(2) = axisPos(2);
            pos(4) = axisPos(4);
            pos(3) = axisPos(4) * pbarRatio;
            pos(1) = (axisPos(3) - pos(3))/2 + axisPos(1);
        end
    end
end

% Convert plot box position to the units used by the axis

hparent = get(h, 'parent');
hfig = ancestor(hparent, 'figure'); % in case in panel or similar
currax = get(hfig, 'currentaxes');

temp = axes('Units', 'Pixels', 'Position', pos, 'Visible', 'off', 'parent', hparent);
set(temp, 'Units', currunit);
pos = get(temp, 'position');
delete(temp);

set(hfig, 'currentaxes', currax);
end




