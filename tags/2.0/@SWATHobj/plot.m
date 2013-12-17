function h = plot(SW,varargin)
% plot instance of SWATHobj in map view
%
% Syntax
%
%     plot(S)
%     plot(S,...)
%     h = ...;
%
% Description
%
%     PLOT creates a map-view plot of different parts of a SWATHobj. The
%     plot may also be color-coded by the z-values, usually elevation, for
%     which either matlab's 'plot' or 'scatter' function can be used. This
%     needs to be specified in the parameters (see below) and affects the
%     speed of the plotting, which can be slow if the SWATHobj contains
%     many data points.
%
% Input arguments
%
%     S     instance of SWATHobj
%
%     Parameter name/value pairs
%
%     'trace'   {true}, false
%     determines if the central trace of the SWATHobj is plotted
%
%     'outline'   {true}, false
%     determines if the outline of the SWATHobj is plotted
%
%     'points'   {true}, false
%     determines if the data points of the SWATHobj are plotted
%
%     'left'   {true}, false
%     determines if the left half (as seen from the direction of the 
%     central line) of the SWATHobj is plotted
%
%     'right'   {true}, false
%     determines if the right half (as seen from the direction of the 
%     central line) of the SWATHobj is plotted
%
%     'legend'   {true}, false
%     determines if a legend is plotted
%
%     'labeldist'   {[]}, scalar, 1 x m vector 
%     adds labels to specific points along the profile. By default, no
%     labels are shown. If a scalar is provided, the profile is divided
%     into segments of equal lengths. If a vector is provided, only the
%     points at the specified distances are labeled. 
%
%     'plotmode'   {'plot'}, 'scatter'
%     switches between Matlab's 'plot' and 'scatter' function for making
%     the plot. The scatter command is uduslly faster when 'colorz' is set
%     to true
%
%     'colorz'   true, {false}
%     optionally colors the data points by the Z value
%
%     'colorrange' 1x2 scalar vector
%     determines the minimum and maximum values of the colormap range.
%     Default values are the minimum and maxmum values of the entire
%     SWATHobj.
%
%     'colormap'   {'jet'}, one of Matlab's builtin colormaps
%     determines the colors used for plotting the data points if 'colorz'
%     is set to true. Note that if plotting on top of, e.g., a DEM produced
%     by imagesc(DEM), the colorbar of the DEM will change accordingly.
%
%     'colormode'   {'normal'}, 'inverse'
%     determines how the colors change with values. Default is 'normal',
%     i.e., low (high) z-values correspond to blue (red) colors in the
%     'jet' colormap.
%
%     'colorbar', {true}, false
%     optionally includes a colorbar; only if 'colorz' is true
%
%     'markersize'   {2}, scalar
%     size of the circular markers used to plot the location of the data 
%     points
%
%
% Output arguments (optional)
%
%     h    handle object of the plotted features. If more than one feature 
%          (trace, outline, points) is plotted, h contains more than one 
%          handle
%
% Examples
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     SW = SWATHobj(DEM,'dx',200,'dy',200);
%     G = gradient8(DEM,'degree');
%     SWG = mapswath(SW,G);
%     figure, plot(SWG,'colorz',true,'plotmode','scatter',...
%         'colorrange',[0 30],'colormap','hot','markersize',10,'legend',false);    
%
%
% Author: Dirk Scherler (scherler[at]caltech.edu)
% Date: June, 2013


cmaps = {'bone','colorcube','cool','copper','flag',...
    'gray','hot','hsv','jet','lines','pink','prism',...
    'spring','summer','white','winter'};


% Parse inputs
p = inputParser;
p.FunctionName = 'plot';
addRequired(p,'SW',@(x) isa(x,'SWATHobj'));
addParamValue(p,'trace',true,@(x) islogical(x))
addParamValue(p,'outline',true,@(x) islogical(x))
addParamValue(p,'points',true,@(x) islogical(x))
addParamValue(p,'left',true,@(x) islogical(x))
addParamValue(p,'right',true,@(x) islogical(x))
addParamValue(p,'legend',true,@(x) islogical(x))
addParamValue(p,'labeldist',[],@(x) isnumeric(x))
addParamValue(p,'plotmode','plot',@(x) ismember(x,{'plot','scatter'}))
addParamValue(p,'colorz',false,@(x) islogical(x))
addParamValue(p,'colormap','jet',@(x) ismember(x,cmaps))
addParamValue(p,'colorrange',[-inf inf],@(x) isnumeric(x))
addParamValue(p,'colormode','normal',@(x) ismember(x,{'normal','inverse'}))
addParamValue(p,'colorbar',true, @(x) islogical(x))
addParamValue(p,'markersize',2,@(x) isnumeric(x))

parse(p,SW,varargin{:});

trace      = p.Results.trace;
outline    = p.Results.outline;
points     = p.Results.points;
left       = p.Results.left;
right      = p.Results.right;
plotlegend = p.Results.legend;
labeldist  = p.Results.labeldist;
plotmode   = p.Results.plotmode;
colorz     = p.Results.colorz;
colmap     = p.Results.colormap;
colrange   = p.Results.colorrange;
colmode    = p.Results.colormode;
colbar     = p.Results.colorbar;
markersize = p.Results.markersize;

if ~left && ~right
    warning('Empty SWATHobj: nothing to plot')
    return;
end

if colorz
    cmap = colormap(colmap)';
    if strcmp(colmode,'inverse'); cmap = flipud(cmap); end
    minz = colrange(1);
    maxz = colrange(2);
    if isinf(minz); minz = min(min(cell2mat(SW.Z))); end
    if isinf(maxz); maxz = max(max(cell2mat(SW.Z))); end
    CX = linspace(minz,maxz,length(cmap));
end

ct =1;
for i = 1 : length(SW.xy0)
    if ~left
        ny = ceil(length(SW.disty{i})/2)+1; % assume all disty are the same
        SW.X{i} = SW.X{i}(ny:end,:);
        SW.Y{i} = SW.Y{i}(ny:end,:);
    elseif ~right
        ny = floor(length(SW.disty{i})/2);
        SW.X{i} = SW.X{i}(1:ny,:);
        SW.Y{i} = SW.Y{i}(1:ny,:);
    end
    
    legct = 1;
    if points
        
        ix = find(~isnan(SW.Z{i}));
        xp = SW.X{i}(ix);
        yp = SW.Y{i}(ix);
        if~isempty(xp)
        if colorz % color-code elevations
            zp = SW.Z{i}(ix);
            % Interpolate colors
            mfcol = [interp1(CX,cmap(1,:),zp),...
                interp1(CX,cmap(2,:),zp),interp1(CX,cmap(3,:),zp)];
            switch plotmode
                case 'plot'
                    for k = 1 : length(zp)
                        ht(ct) = plot(xp(k),yp(k),'o','color',mfcol(k,:),'Markersize',markersize); hold on
                    end
                case 'scatter'
                    ht(ct) = scatter(xp,yp,markersize,mfcol);
            end
        else
            mfcol = [0 0 0];
            switch plotmode
                case 'plot'
                    ht(ct) = plot(xp,yp,'o','color',mfcol,'Markersize',markersize); hold on
                case 'scatter'
                    ht(ct) = scatter(xp,yp,markersize,mfcol);
            end
        end
        end
        legstr(legct) = {'Points'};
        ct=ct+1; legct = legct+1;
    end
    if trace
        ht(ct) = plot(SW.xy0{i}(:,1),SW.xy0{i}(:,2),'go','Markersize',2); hold on
        legstr(legct) = {'Original trace'}; 
        ct=ct+1; legct = legct+1;
        ht(ct) = plot(SW.xy{i}(:,1),SW.xy{i}(:,2),'g-','Markersize',2);
        legstr(legct) = {'Resamp./interp. trace'};
        ct=ct+1; legct = legct+1;
    end
    if outline
        IM = ~isnan(SW.Z{i});
        B = bwboundaries(IM,4);
        for k = 1 : length(B)
            ix = sub2ind(size(IM),B{k}(:,1),B{k}(:,2));
            x_outline = SW.X{i}(ix);
            y_outline = SW.Y{i}(ix);
            ht(ct) = plot(x_outline,y_outline,'r-');
        end
        legstr(legct) = {'Outline'};
        ct=ct+1; legct = legct+1;
    end
    
    if ~isempty(labeldist)
        xp = SW.xy{i}(:,1);
        yp = SW.xy{i}(:,2);
        dp = SW.distx{i};
        if isscalar(labeldist)
            lp = (0:labeldist:max(dp))';
        elseif isvector(labeldist)
            lp = labeldist;
        end
        % interpolate label positions
        lx = interp1(dp,xp,lp);
        ly = interp1(dp,yp,lp);
        % add labels
        textstr = cellstr(num2str(lp));
        hold on, plot(lx,ly,'wv')
        text(lx+2*SW.dx,ly,textstr,'color',[1 1 1]); % offset=2xSW.dx
    end
end

drawnow
axis equal
if plotlegend
    legend(ht(1:length(legstr)),legstr);
end

if colorz && colbar
    hc = colorbar;
    yt = get(hc,'YTick');
    set(hc,'YTickLabel',(yt.*(maxz-minz)+minz)');
end


if nargout>0
    h = ht;
end

end