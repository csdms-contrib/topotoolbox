function slopeareatool(FD,DEM,varargin)

%SLOPEAREATOOL Interactively create slope area plots and fit power laws
%
% Syntax
%     
%     slopeareatool(FD,DEM)
%     slopeareatool(FD,DEM,'pn','pv',...)
%
% Description
%
%     SLOPEAREATOOL is an interactive tool to map streams and flowpaths and
%     simultaneously plot area vs. gradient. Mapping can be done using up
%     to five different plot colors so that user-defined groups of streams
%     can be created. A curve-fitting scheme allows to fit group-specific
%     power laws (S = beta(1) A^(beta(2)) and thus enables to compare these
%     groups based on the derived parameters beta.
%
% Input arguments
%
%     FD    Flow direction (FLOWobj)
%     DEM   digital elevation  model (GRIDobj)
%
%     Parameter name/value pairs
%  
%     minarea        minimum upslope area for stream grid derivation in 
%                    map units^2. If set, channelhead are snapped to the
%                    streams.
%     maxarea        maximum upslope area in map units^2. Must not be 
%                    less than minarea (may be set if only interested in
%                    tributary streams)
%     nrbins         nr of bins ({500}) between minium and maximum of 
%                    upslope grid. If minarea/maxarea are set, the bins are
%                    distributed with logarithmic spacing between minarea 
%                    and maxarea.
%     mingradient    minimum gradient. Since negative or zero gradients
%                    will not be displayed in loglog space, they will be
%                    increased to minimum gradient ({0.0001})
%     fitmethod      see help slopearea for description
%
%
% See also: slopearea, chiplot
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 19. June, 2013


p = inputParser;
p.FunctionName = 'slopeareatool';

addRequired(p,'FD',@(x) isa(x,'FLOWobj'));
addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));

validfitmethods  = {'ls','lad'};
addParamValue(p,'minarea',0,@(x) isscalar(x));
addParamValue(p,'maxarea',inf,@(x) isscalar(x));
addParamValue(p,'nrbins',50,@(x) isscalar(x));
addParamValue(p,'mingradient',0.0001,@(x) isscalar(x) && x>0);
addParamValue(p,'fitmethod','ls');

parse(p,FD,DEM,varargin{:});

% required
FD         = p.Results.FD;
DEM        = p.Results.DEM;
minarea    = p.Results.minarea;
maxarea    = p.Results.maxarea;
nrbins     = p.Results.nrbins;

fitmethod = validatestring(p.Results.fitmethod,validfitmethods);

% validate alignment of FD and DEM
validatealignment(FD,DEM)

% If you set fastindexing to true then downstream processing is much faster
% for tools that do computations along single flow paths such as
% flowpathextract
FD.fastindexing = true;

% create figure

scrsz = get(0,'ScreenSize');
% pos = [left, bottom, width, height]
hFig = figure('OuterPosition',[1/8*scrsz(3) 1/3*scrsz(4) 2/4*scrsz(3) 2/3*scrsz(4)],...
              'MenuBar','none',...
              'NumberTitle','off',...
              'Name','Main');
          
          
% create menu
hMenuView    = uimenu(hFig,'Label','Change color'); 
hButtonMenuColor = uimenu(hMenuView,'Label','Marker Color'); 
colors = 'kbgyr';
colornames = {'black' 'blue' 'green' 'yellow' 'red'};
% Default line color
props.linecolor = 'k';

hButtonColors(1) = uimenu(hButtonMenuColor,'Label',colornames{1},'Checked','on','Callback',@(src,event) changelinecolor(src,event)); 
hButtonColors(2) = uimenu(hButtonMenuColor,'Label',colornames{2},'Checked','off','Callback',@(src,event) changelinecolor(src,event));
hButtonColors(3) = uimenu(hButtonMenuColor,'Label',colornames{3},'Checked','off','Callback',@(src,event) changelinecolor(src,event));
hButtonColors(4) = uimenu(hButtonMenuColor,'Label',colornames{4},'Checked','off','Callback',@(src,event) changelinecolor(src,event));
hButtonColors(5) = uimenu(hButtonMenuColor,'Label',colornames{5},'Checked','off','Callback',@(src,event) changelinecolor(src,event));
hButtonClear = uimenu(hMenuView,'Label','Clear','Callback',@(src,event) clearvectorplots);  

hMenuFit =   uimenu(hFig,'Label','Curve fitting','Enable','off');

hMenuExport  = uimenu(hFig,'Label','Export','Enable','off');
      
          
          
% calculate hillshade
RGB  = imageschs(DEM,DEM);

% calculate terrain attributes
% Upslope area (we will use pixel units and switch to map units when plotting)
A = flowacc(FD);
% gradient
DEM = imposemin(FD,DEM,0.0000001);
G = FLOWobj2gradient(FD,DEM);

% Set minarea and maxarea to pixel units
minarea = minarea/(A.cellsize^2);
maxarea = maxarea/(A.cellsize^2);


%% Maximum upslope area set?
if isinf(maxarea);
    maxarea = max(A);
elseif maxarea < minarea
    error('TopoToolbox:slopeareatool','maxarea must be larger than minarea')
else
    FD.ixcix(A.Z>maxarea) = 0;  
end    


%% Snapping
if minarea > 0;
    % calculate stream grid
    if minarea < 10;
        warning('TopoToolbox:slopeareatool','minimum area is quite small')
    end
    W = A.Z>(minarea) & A.Z<(maxarea);
    
    RGB(repmat(W,[1 1 3])) = 150;    
    [~,SNAPRASTER] = bwdist(W,'q');

    snap = true;
else
    minarea = 1;
    snap = false;
end
    
%% Plot hillshade and magnifying tool   
% create empty axes in figure
% hAx = axes('parent',hFig);
hAx = imgca(hFig);

% show hillshade using imshow
hIm = image(RGB,'parent',hAx);
% create an instance of imscrollpanel and use the API
hPanel = imscrollpanel(hFig,hIm);
api    = iptgetapi(hPanel);
% set to initial magnification
mag    = api.findFitMag();
api.setMagnification(mag);

% create magnification box in another window
hMagBox = immagbox(hFig,hIm);
pos     = get(hMagBox,'Position');
set(hMagBox,'Position',[0 0 pos(3) pos(4)])
hFigOV = imoverview(hIm);

% create slope area figure
hAxProfiles = createSlopeAreaFigure;


%% distribute figures on screen
% distFig('Position','L','Only',[hFig],'Rows',1)
% distFig('Position','R','Only',[hFigOV hFigProfiles],'Rows',2)

%% add callbacks
enterFcn = @(figHandle, currentPoint)...
       set(hFig, 'Pointer','crosshair');
iptSetPointerBehavior(hIm,enterFcn);
iptPointerManager(hFig);

xlim = get(hIm,'XData');
ylim = get(hIm,'YData');
X = 1:xlim(2);
Y = 1:ylim(2);

set(hFig,'WindowButtonDownFcn',@PressLeftButton);

%% Predefine variables
IXchannelhead = 0;
counter = 0;
color   = [];
IXchannel = {};
beta  = zeros(2,numel(colors));
hPlot = [];
hPlotProfiles = [];
hPlotFit = inf(numel(colors),1);
hLegendFit = [];
hButtonFit = [];
hButtonExport = [];

%% Create bins
% binning
aedges  = logspace(floor(log10(minarea-.01)),ceil(log10(maxarea+1)),nrbins+1)';
acenters = aedges + [diff(aedges)/2;0];




    %% Subfunctions
    function PressLeftButton(src,evt)
        
        % get the axis position
        cp = get(hAx,'CurrentPoint');
        
        % get pixel locations
        cp = cp(1,1:2);
        pixelx = round(axes2pix(xlim(2), xlim, cp(1)));
        pixely = round(axes2pix(ylim(2), ylim, cp(2)));
        
        % check if current point is located inside axis
        if pixelx < xlim(1) || pixelx > xlim(2) ...
                || pixely < ylim(1) || pixely > ylim(2)
            
        else
            % call function
            IXchannelhead(:) = sub2ind(DEM.size,pixely,pixelx);
            counter = counter + 1;
            
            % snap to stream
            if snap
                IXchannelhead(:) = SNAPRASTER(IXchannelhead);
            end
            
            PlotProfiles
        end
        
    end



    function PlotProfiles
        % get coordinates of mouse click
        if A.Z(IXchannelhead) > maxarea;
            counter = counter-1;
            return
        elseif FD.ixcix(IXchannelhead) == 0
            counter = counter -1;
            return
        end
        
        % get indices of flow path
        IXchannel{counter} = flowpathextract(FD,IXchannelhead);
        FD.ixcix(IXchannel{counter}) = 0;
        
        % color
        color(counter) = strfind(colors,props.linecolor);
        
        if counter == 1;
            set(hMenuFit,'Enable','on')
            hButtonFit = uimenu(hMenuFit,'Label',['Fit to ' colornames{color(counter)} ' data'],'Callback',@(src,event) fitline(src,event,color(counter)));
            
            set(hMenuExport,'Enable','on')
            hButtonExport = uimenu(hMenuExport,'Label',['Export ' colornames{color(counter)} ' streams as STREAMobj to workspace'],'Callback',@(src,event) exporttoworkspace(src,event,color(counter))); 
        elseif ~ismember(color(end),color(1:end-1));
            hButtonFit(end+1) = uimenu(hMenuFit,'Label',['Fit to ' colornames{color(counter)} ' data'],'Callback',@(src,event) fitline(src,event,color(counter)));
            hButtonExport(end+1) = uimenu(hMenuExport,'Label',['Export ' colornames{color(counter)} ' streams as STREAMobj to workspace'],'Callback',@(src,event) exporttoworkspace(src,event,color(counter)));
        end
        
        % bin area values
        [~,ix]  = histc(A.Z(IXchannel{counter}),aedges);
        ix = ix(:);
        I  = ix>0;
        I(end) = false;
        IXchannel{counter}(~I) = [];
        
        try
            g  = accumarray(ix(I),double(G.Z(IXchannel{counter})),[numel(aedges) 1],@(x) mean(double(x)),nan);
        catch %#ok<CTCH>
            counter = max(counter - 1,0);
            return
        end
        
        
        % plot flow path
        hold(hAx,'on')
        [r,c] = ind2sub(DEM.size,IXchannel{counter});
        % let the new flow path blink
        hPlot(counter) = plot(hAx,X(c),Y(r),props.linecolor,'LineWidth',2);
        pause(.1)
        set(hPlot(counter),'Color','w');
        pause(.2)
        set(hPlot(counter),'Color',props.linecolor);
        drawnow
        hold(hAx,'off');
        
        % Slope area plot
        if ~ishandle(hAxProfiles);
            hAxProfiles = createSlopeAreaFigure;
        end
        
        hold(hAxProfiles,'on')
        hPlotProfiles(counter) = loglog(hAxProfiles,acenters*(A.cellsize^2),g,'s',...
            'MarkerSize',8,...
            'MarkerEdgeColor',props.linecolor,...
            'MarkerFaceColor','none');%props.linecolor
        hold(hAxProfiles,'off');
        drawnow
        
        
    end

    function changelinecolor(src,event)
        
        I = src == hButtonColors;
        
        props.linecolor = colors(I);
        
        h = findobj(hButtonColors,'Checked','on');
        set(h,'Checked','off');
        set(src,'Checked','on');
        
    end

    function fitline(src,event,colornum)
        
        switch get(src,'Checked')
            case 'off'
                
                linecolor = colors(colornum);
                
                % fit line
                IXchannel = IXchannel(:);
                IX = cell2mat(IXchannel(color == colornum));
                IX = unique(IX);
                
                S  = convert2STREAMobj(FD,IX);
                SA = slopearea(S,DEM,A,...
                    'plot',false,...
                    'mingradient',p.Results.mingradient,...
                    'fitmethod',fitmethod,...
                    'streamgradient','forward',...
                    'gradaggfun','mean');
                
                beta(:,colornum) = [SA.ks SA.theta]'; 
                
                geval = beta(1,colornum)*((acenters*(A.cellsize^2)).^beta(2,colornum));
                
               
                hold(hAxProfiles,'on')
                hPlotFit(colornum) = plot(hAxProfiles ,acenters*(A.cellsize^2),geval,linecolor,'LineWidth',1.5);
                hold(hAxProfiles,'off')
                
               
                
                % redraw legend         
                if ~isempty(hLegendFit)
                else
                    delete(hLegendFit);
                    hLegendFit = [];
                end
                    
                legplothandle = [];
                M = {};
                for r = 1:numel(colors)
                    if ishandle(hPlotFit(r))
                        legplothandle(end+1) = hPlotFit(r);
                        M{end+1} = ['S = ' num2str(beta(1,r)) ' A^{' num2str(beta(2,r)) '}'];
                    end
                end
                if any(legplothandle)
                    hLegendFit = legend(legplothandle,M,'location','northeast');
                end
                
                
                set(src,'Checked','on');
                
                
            case 'on'
                set(src,'Checked','off');
                
                delete(hPlotFit(colornum))
                hPlotFit(colornum) = inf;
                
                % redraw legend         
                if ~isempty(hLegendFit)
                else
                    delete(hLegendFit);
                    hLegendFit = [];
                end
                    
                legplothandle = [];
                M = {};
                for r = 1:numel(colors)
                    if ishandle(hPlotFit(r))
                        legplothandle(end+1) = hPlotFit(r);
                        M{end+1} = ['S = ' num2str(beta(1,r)) ' A^{' num2str(beta(2,r)) '}'];
                    end
                end
                if any(legplothandle)
                    hLegendFit = legend(legplothandle,M,'location','northeast');
                end
                
        end
    end

    function exporttoworkspace(src,event,colornum)
        linecolor = colors(colornum);
        % fit line
        IXchannel = IXchannel(:);
        IX = cell2mat(IXchannel(color == colornum));
        IX = unique(IX);
        
        S  = convert2STREAMobj(FD,IX);
        
        prompt = {'Enter variable name:'};
        title = 'Export';
        lines = 1;
        def = {['S' linecolor]};
        answer = inputdlg(prompt, title, lines, def);
        if ~isempty(answer) && isvarname(answer{1})
            assignin('base',answer{1},S);
        else
            return
        end
    end
        
        

    function clearvectorplots
        delete(hPlot);
        hPlot = [];
        if ishandle(hPlotProfiles)
            delete(hPlotProfiles)
        end
        hPlotProfiles = [];
        delete(hPlotFit(ishandle(hPlotFit)))
        hPlotFit = inf(numel(colors),1);
        counter = 0;
        color   = [];
        FD.fastindexing = true;
        
        if snap
            FD.ixcix(~W) = 0;
        end
        
        IXchannel = {};
        
        set(hMenuFit,'Enable','off');
        delete(hButtonFit);
        hButtonFit = [];
        
        set(hMenuExport,'Enable','off');
        delete(hButtonExport);
        hButtonExport = [];
        
        if ~isempty(hLegendFit)
            delete(hLegendFit);
            hLegendFit = [];
        end
    end


    function S = convert2STREAMobj(FD,IX)
        
        WW = DEM;
        WW.Z = false(DEM.size);
        WW.Z(IX) = true;
        
        S  = STREAMobj(FD,WW);
    end
        
    function hAxProfiles = createSlopeAreaFigure
        %% create figure for profiles
        hFigProfiles = figure('name','Slope-Area Plot'); 
        hAxProfiles = axes('Parent',hFigProfiles,'Xscale','log','Yscale','log','box','on');
        xlabel(hAxProfiles,'Area [m^2]')
        ylabel(hAxProfiles,'Gradient [mm^{-1}]')
    end
        
        
        

end
        

