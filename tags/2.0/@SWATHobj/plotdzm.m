function h = plotdzm(SW,M,varargin)
% plot instance of SWATHobj in profile view, color-coded based on GRIDobj
%
% Syntax
%
%     plotdzm(SW,M)
%     plotdzm(SW,M,...)
%     h = ...;
%
% Description
%
%     PLOTDZM creates a profile-view plot of a SWATHobj, showing the
%     z-values as a function of distance along the profile of the SWATHobj
%     and color-coded by the z-values in the GRIDobj M.
%
%     For display the function can use either Matlab's 'plot' or 'scatter' 
%     function or create an image. This needs to be specified in the 
%     parameters (see below) and affects the speed of the plotting, which 
%     can be slow if the SWATHobj contains many data points. In this case,
%     creating an image is fastest.
%
% Input arguments
%
%     S     instance of SWATHobj
%     M     instance of GRIDobj
%
%     Parameter name/value pairs
%
%     'ix'   {[]}, 1 x m vector
%     list of indices to the SWATHobj that specifies which entries to plot.
%
%     'left'   {true}, false
%     determines if the left half (as seen when looking along the direction
%     of the central line) of the SWATHobj is plotted.
%
%     'right'   {true}, false
%     determines if the right half (as seen when looking along the 
%     direction of the central line) of the SWATHobj is plotted.
%
%     'dist'   {'x'}, 'y'
%     plot data either along the x-direction (default) of the SWATHobj or 
%     along the y-direction (across the SWATHobj).
%
%     'distadjust'   {0}, scalar
%     allows shifting the x-axis by a scalar value. This is useful when
%     alligning SWATHobj's that were obtained along, e.g., a certain reach of
%     a drainage network, with the STREAMobj that it was derived from.
%
%     'merge'    true, {false}
%     determines if different components of a SWATHobj will be plotted in
%     one or separate figures.
%
%     'zmode'   {'absolute'}, 'relative'
%     assigns elevations along the SWATHobj either in absolute values
%     (default), or relative to the minimum value across each line in the
%     SWATHobj (relative). Note that this can be inaccurate if the points
%     across the SWATHobj miss to sample the lowest elevation present.
%
%     'sortm'   {'descend'}, 'ascend', 'none'
%     determines if and how the values within each component of a SWATHobj
%     are sorted before being plotted. This parameter affects the
%     visibility of specific populations of the data points if the markers
%     overlap in the plot.
%
%     'colormap'   {'jet'}, one of Matlab's builtin colormaps
%     determines the colors used for plotting the data points if 'colorz'
%     is set to true. Note that if plotting on top of, e.g., a DEM produced
%     by imagesc(DEM), the colorbar of the DEM will change accoringly.
%
%     'colormode'   {'normal'}, 'inverse'
%     determines how the colors change with values. Default is 'normal',
%     i.e., low (high) m-values would correspond to blue (red) colors in 
%     the 'jet' colormap.
%
%     'colorrange'   1x2 vector
%     determines the minimum and maximum values of the colormap range.
%     Default values are the minimum and maxmum values of the entire
%     SWATHobj.
%
%     'colorbar'   {true}, false
%     adds a colorbar to the plot.
%
%     'plotmode'   {'plot'}, 'scatter', 'image'
%     defines how the data points are plotted. If the number of data is
%     rather small 'plot' and 'scatter' may produce a figure that can be
%     exported as editable vector graphics. This feature breaks down at
%     larger numbers of data points and the time required for building the
%     figure can get very long. In such cases it may be useful to 
%     interpolate an 'image'.
%
%     'markersize'   {2}, scalar
%     size of the circular markers used to plot the location of the data 
%     points. Only applicable if 'plotmode' set to 'plot' or 'scatter'.
%
%     'xstretch' and 'ystretch'    {1 and 1}, scalar
%     if 'plotmode' is set to 'image', the x- and y-values can be stretched
%     by the corresponding factors to increase the resolution of the image.
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
%     FD = FLOWobj(DEM,'preprocess','carve');
%     FA = flowacc(FD);
%     S = STREAMobj(FD,FA>1e6/FA.cellsize/FA.cellsize);
%     S = klargestconncomps(S);
%     S = trunk(S);
%     SW = SWATHobj(DEM,S,'dx',100,'dy',100,'width',2000,'smooth',11,'plot',false);
%     imagesc(DEM), hold on, plot(SW)
%     G = gradient8(DEM,'degree');
%     plotdzm(SW,G,'colorrange',[0 35]); 
%
%
% Author: Dirk Scherler (scherler[at]caltech.edu)
% Date: June, 2013


cmaps = {'bone','colorcube','cool','copper','flag',...
    'gray','hot','hsv','jet','lines','pink','prism',...
    'spring','summer','white','winter'};

% Parse inputs
p = inputParser;
p.FunctionName = 'plotdzm';
addRequired(p,'SW',@(x) isa(x,'SWATHobj'));
addRequired(p,'M',@(x) isa(x,'GRIDobj'));
addParamValue(p,'ix',[],@(x) isvector(x))
addParamValue(p,'left',true,@(x) islogical(x))
addParamValue(p,'right',true,@(x) islogical(x))
addParamValue(p,'dist','x',@(x) ismember(x,{'x','y'}))
addParamValue(p,'distadjust',0,@(x) isnumeric(x))
addParamValue(p,'merge',false,@(x) islogical(x))
addParamValue(p,'zmode','absolute',@(x) ismember(x,{'absolute','relative'}))
addParamValue(p,'colormap','jet',@(x) ismember(x,cmaps))
addParamValue(p,'colormode','normal',@(x) ismember(x,{'normal','inverse'}))
addParamValue(p,'colorrange',[-inf,inf],@(x) isnumeric(x))
addParamValue(p,'colorbar',true, @(x) islogical(x))
addParamValue(p,'sortm','descend',@(x) ismember(x,{'ascend','descend','none'}))
addParamValue(p,'plotmode','scatter',@(x) ismember(x,{'plot','scatter','image'}))
addParamValue(p,'xstretch',1, @(x) isnumeric(x) & x>0)
addParamValue(p,'ystretch',1, @(x) isnumeric(x) & x>0)
addParamValue(p,'markersize',4,@(x) isnumeric(x))


parse(p,SW,M,varargin{:});

% parameters
thisix     = p.Results.ix;
left       = p.Results.left;
right      = p.Results.right;
dist       = p.Results.dist;
distadjust = p.Results.distadjust;
merge      = p.Results.merge;
zmode      = p.Results.zmode;
colmap     = p.Results.colormap;
colmode    = p.Results.colormode;
colrange   = p.Results.colorrange;
colbar     = p.Results.colorbar;
sortm      = p.Results.sortm;
plotmode   = p.Results.plotmode;
markersize = p.Results.markersize;
xstretch   = p.Results.xstretch;
ystretch   = p.Results.ystretch;


warning('off')
[SWG] = mapswath(SW,M);
warning('on')
    
if ~left && ~right
    warning('Empty SWATHobj: nothing to plot')
    return;
end

minm = colrange(1);
maxm = colrange(2);
if isinf(minm); minm = min(min(cell2mat(SWG.Z))); end
if isinf(maxm); maxm = max(max(cell2mat(SWG.Z))); end

if ~strcmp(plotmode,'image')
    hf = figure('visible','off');
    cmap = colormap(colmap)';
    % cb = colorbar;
    % cdata = get(get(cb,'Children'),'Cdata');
    close(hf)
    if strcmp(colmode,'inverse'); cmap = fliplr(cmap); end
    CX = linspace(minm,maxm,length(cmap));
end


if isempty(thisix)
    thisix = 1 : length(SW.xy0);
end

hout = zeros(size(thisix));
for m = 1 : length(thisix)
    
    i = thisix(m);
    
    if ~merge
        figure;
    end
    
    hout(m) = gcf;
    
    if ~left
        ny = ceil(length(SW.disty{i})/2)+1;
        SW.Z{i} = SW.Z{i}(ny:end,:);
        SWG.Z{i} = SWG.Z{i}(ny:end,:);
    elseif ~right
        ny = floor(length(SW.disty{i})/2);
        SW.Z{i} = SW.Z{i}(1:ny,:);
        SWG.Z{i} = SWG.Z{i}(1:ny,:);
    end
    
    if strcmp(zmode,'relative')
        minm = min(SW.Z{i},[],1);
        SW.Z{i} = SW.Z{i}-repmat(minm,length(SW.disty{i}),1);
    end

    [zr,zc] = size(SW.Z{i});
    switch dist
        case 'x'
            d = repmat(SW.distx{i}',zr,1);
            swdx = SW.dx;
        case 'y'
            d = repmat(SW.disty{i},1,zc);
            swdx = SW.dy;
    end
    z = reshape(SW.Z{i},numel(SW.Z{i}),1);
    s = reshape(SWG.Z{i},numel(SWG.Z{i}),1);
    d = reshape(d,numel(d),1);
    
    d = d(~isnan(z));
    s = s(~isnan(z));
    z = z(~isnan(z));
    
    
    % Sort by slope
    if ismember(sortm,{'ascend','descend'})
        [s,ix] = sort(s,sortm);
        z = z(ix);
        d = d(ix);
    end
        
    % Interpolate colors
    if ~strcmp(plotmode,'image')
        mfcol = [interp1(CX,cmap(1,:),s),...
            interp1(CX,cmap(2,:),s),interp1(CX,cmap(3,:),s)];
    end
    
    switch plotmode
        case 'plot'
            for k = 1 : length(z)
                plot(d(k)+distadjust,z(k),'s','color',mfcol(k,:),'Markerfacecolor',mfcol(k,:),'Markersize',markersize), hold on
            end
            
        case 'scatter'
            scatter(d+distadjust,z,markersize,mfcol,'square','filled'), hold on
            
        case 'image'
            x = d+distadjust;
            y = z;
            v = s;
            
            % apply colorrange thresholds
            x = x(v>=minm);
            y = y(v>=minm);
            v = v(v>=minm);
            v = min(v,maxm);
            y = round(y.*ystretch); % stretch elevation-values by 10
            x = round(x.*xstretch); % stretch elevation-values by 10
            
            F = TriScatteredInterp(x,y,v);
            
            xy = unique([x,y],'rows');
            qv = F(xy(:,1),xy(:,2));
            qx = floor(min(x)):swdx:ceil(max(x));
            qy = floor(min(y)):ceil(max(y));
            IM = GRIDobj(qx./swdx,qy,zeros(length(qy),length(qx)));
            ix = coord2ind(IM,xy(:,1)./swdx,xy(:,2));
            IM.Z(ix) = qv;
            
%             SE = strel('square',3)./9;
%             IM.Z = imdilate(IM.Z,SE);
%             SE = ones(5)./25;
%             IM.Z = conv2(IM.Z,SE);
            
            hout = imagesc(qx./xstretch,flipud(qy'./ystretch),IM.Z,[minm maxm]);
            set(gca,'YDir','normal');
            hc = colormap(colmap);
            if strcmp(colmode,'inverse')
                colormap(flipud(hc));
            end
            axis normal
    end
    
end % loop over SWATHobj entries
drawnow
xlabel(sprintf('Distance along profile (%s)',SW.xyunit))
ylabel(sprintf('Z (%s)',SW.zunit))

if colbar && ~strcmp(plotmode,'image')
    colormap(colmap);
    hc = colorbar;
    if strcmp(colmode,'inverse');
        set(get(hc,'Children'),'CData',(64:-1:1)');
    end
    yt = get(hc,'YTick');
    set(hc,'YTickLabel',(yt.*(maxm-minm)+minm)');
elseif colbar && strcmp(plotmode,'image')
    colorbar;
end

if nargout==1;
    h = unique(hout);
end

