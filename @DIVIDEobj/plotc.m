function varargout = plotc(D,val,varargin)
%PLOTC   Create colored plot of drainage divide network
%
% Syntax
%
%     plotc(D,GRID)
%     plotc(D,val)
%     plotc(D,...,'pn','pv')
%     h = plot(D)
%
% Description
%
%     PLOTC creates a colored plot of the divide network. The divide network
%     is shown as a collection of points, one for each divide edge. The
%     points can vary in size by the along-divide distance or by the divide
%     order; they can vary in color by attributes provided either in form
%     of a GRIDobj or as a list of values, with the same number of elements 
%     as divide edges. The appearance of the colors can also be modified
%     by using the HSV color model, in which case the saturation and
%     'value' are linked to another GRIDobj or list of values.
%     
% Input arguments
%
%     D      Drainage divide network (class DIVIDEobj)
%     GRID   GRIDobj of the same extent as the DEM, from which the divides
%            were obtained.
%     val    vector the same size as the field "D.IX" (see function
%            getvalue)
%
%     Parameter name/value pairs
%
%     'divsize'    size of the divides {'distance'(default),'order','none'}
%     'size'       minimum marker size (default = [2 50])
%     'valuefct'   function to compute divide color value from adjacent
%                  pixels if GRIDobj is provided {'mean'(default),'min',
%                  'max','diff','normdiff'}
%     'caxis'      colorrange limits
%     'divlim'     limit plot of divides {'distance'(default),'order','none'}
%     'limit'      lower limit (default = [-inf inf])
%     'hsvdata'    data to set saturation and value in hsv color space
%     'hsvrange'   minimum saturation and value (default = [0 1])
%     'hsvfct'     function to compute divide color value from adjacent
%                  pixels if GRIDobj is provided {'mean'(default),'min',
%                  'max','diff','normdiff'}
%
% Output arguments
%
%     h  graphics handle
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','c');
%     ST = STREAMobj(FD,flowacc(FD)>1000);
%     D = DIVIDEobj(FD,ST);
%     D = divorder(D,'topo');
%     plotc(D,D.distance./1e3,'limit',[1000 inf])
%     box on
%     axis image
%     hc = colorbar;
%     hc.Label.String = 'Divide distance, d_d (km)';
%     
%
% See also: DIVIDEobj, DIVIDEobj/plot, getvalue
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: April 2020


%% Parse inputs
p = inputParser;
p.KeepUnmatched = true;
addParameter(p,'divsize','distance',...
    @(x) ischar(validatestring(x,{'distance','order','none'})));
addParameter(p,'size',[2 50],@(x) isnumeric(x));
addParameter(p,'valuefct','mean',...
    @(x) ischar(validatestring(x,{'min','max','mean','diff','normdiff'})));
addParameter(p,'caxis',[]);
addParameter(p,'divlimit','distance',...
    @(x) ischar(validatestring(x,{'distance','order','none'})));
addParameter(p,'limit',[-inf inf],@(x) isnumeric(x));
addParameter(p,'hsvdata',[],@(x) isnumeric(x));
addParameter(p,'hsv',[0 1],@(x) isnumeric(x));
addParameter(p,'hsvfct','mean',...
    @(x) ischar(validatestring(x,{'min','max','mean','diff','normdiff'})));

% Parse
parse(p,varargin{:});
hsvdata = p.Results.hsvdata;

siz = p.Results.size;

if strcmp(p.Results.divsize,'none')
    siz(:) = min(siz);
end

[x,y] = ind2coord(D,D.IX);

if isa(val,'GRIDobj')
    [q,r] = getvalue(D,val);
    
    switch p.Results.valuefct
        case 'mean'
            val = nanmean([q,r],2);
        case 'max'
            val = nanmax([q,r],[],2);
        case 'min'
            val = nanmin([q,r],[],2);
        case 'diff'
            val = abs(diff([q,r],1,2));
        case 'normdiff'
            val = abs(diff([q,r],1,2))./nanmean([q,r],2);
    end
end


if not(isempty(hsvdata)) 
    if isa(hsvdata,'GRIDobj')
        [q,r] = getvalue(D,hsvdata);

        switch p.Results.hsvfct
            case 'mean'
                hsvdata = nanmean([q,r],2);
            case 'max'
                hsvdata = nanmax([q,r],[],2);
            case 'min'
                hsvdata = nanmin([q,r],[],2);
            case 'diff'
                hsvdata = abs(diff([q,r],1,2));
            case 'normdiff'
                hsvdata = abs(diff([q,r],1,2))./nanmean([q,r],2);
        end
    end    
end


if not(isempty(p.Results.caxis))
    ix = val>max(p.Results.caxis);
    val(ix) = max(p.Results.caxis);
    ix = val<min(p.Results.caxis);
    val(ix) = min(p.Results.caxis);
end


if strcmp(p.Results.divsize,'order') && not(isempty(D.order))
    divsiz = D.order;
elseif strcmp(p.Results.divsize,'distance') && not(isempty(D.distance))
    divsiz = D.distance;
else
    divsiz = ones(size(val));
end

if strcmp(p.Results.divlimit,'order')
    ix = D.order>=min(p.Results.limit) & ...
        D.order<=max(p.Results.limit);
elseif strcmp(p.Results.divlimit,'distance')
    ix = D.distance>=min(p.Results.limit) & ...
        D.distance<=max(p.Results.limit);
else
    ix = true(size(val));
end

% scale size
a = range(siz)/range(divsiz(ix));  % ISSUE: can be NaN!
b = siz(1)-min(divsiz(ix)).*a;
ms = a.*divsiz+b;

if isempty(hsvdata)
    h = scatter(x(ix),y(ix),ms(ix),val(ix),'filled');
else
    % use hsv colormap
	cmap = hsv(256);
    cmap = rgb2hsv(cmap);
    % stretch colormap according to minimum and maximum values
    minval = min(val(ix));
    maxval = max(val(ix));
    X = linspace(minval,maxval,length(cmap));
    % interpolate hue
    col = zeros(length(val),3);
	col(:,1) = interp1(X,cmap(:,1),val);
    col(:,2) = interp1(X,cmap(:,2),val);
    col(:,3) = interp1(X,cmap(:,3),val);
    % use hsvdata to set saturation and value of hsv colors
    hsvdata = hsvdata-min(hsvdata);
    hsvdata = hsvdata./max(hsvdata);
    hsvdata = hsvdata-min(p.Results.hsv);
    hsvdata(hsvdata<0) = 0;
    hsvdata = hsvdata./max(hsvdata);
    hsvdata = hsvdata./max(p.Results.hsv);
    hsvdata(hsvdata>1) = 1; 
    
    col(:,2) = hsvdata;
    col(:,3) = hsvdata;
    col = hsv2rgb(col);
    h = scatter(x(ix),y(ix),ms(ix),col(ix,:),'filled');
    colormap hsv;
end

if not(isempty(p.Results.caxis))
    hca = gca;
    hca.CLim = p.Results.caxis;
end

if nargout==1
    varargout{1} = h;
end




end

