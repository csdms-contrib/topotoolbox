function varargout = getvalue(D,GRID,varargin)
%GETVALUE    get grid values adjacent to divides
%
% Syntax
% 
%     [p,q] = getvalue(D,GRID)
%     [r] = getvalue(D,GRID)
%     [r] = getvalue(D,GRID,fct)
%     [p,q,pix,qix] = getvalue(D,GRID);
%
% Description
%
%     GETVALUE provides the GRIDobj values of pixels that are adjacent to
%     divide edges. 
%     
%
% Input arguments
%
%     D       drainage divide network (class DIVIDEobj)
%     GRID    grid of the same area as D (class GRIDobj)
%     fct     name of function used to aggregate the values in p,q to
%             compute r {'min','max','mean'(default),'diff','normdiff'}
% 
% Output arguments
%     
%     p,q     GRIDobj values of pixels adjacent to divides
%     r       GRIDobj values of pixels adjacent to divides, aggregated into
%             one value per divide edge
% 
% Optional output
%
%    pix,qix  linear indices of the p,q data in the GRIDobj
%    
%
% Example
%     
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','c');
%     ST = STREAMobj(FD,flowacc(FD)>1000);
%     D = DIVIDEobj(FD,ST);
%     D = divorder(D,'topo');
%     D = divdist(D);
%     DZ = vertdistance2stream(FD,ST,DEM);
%     DZ.Z(isinf(DZ.Z)) = nan;
%     [p,q] = getvalue(D,DZ);
%     dz = abs(diff([p,q],1,2));
%     plotc(D,dz,'caxis',[0 300],'limit',[1000 inf])
%     colorbar
%     title('Across-divide difference in hillslope relief')
%     
%     
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: December 2018


% Get vectors 
[x,y] = ind2coord(D,vertcat(D.IX));

% Find pixels on either side of ridgeline
x1 = [NaN; x(1:end-1)];
x2 = [NaN; x(2:end)];
y1 = [NaN; y(1:end-1)];
y2 = [NaN; y(2:end)];
dx = x1-x2;
dy = y1-y2;
hcs = D.cellsize/2;
ix = dx==0; % vertical link
iy = dy==0; % horizontal link
meanx = (x1+x2)./2;
meany = (y1+y2)./2;
px = meanx + hcs.*ix;
qx = meanx - hcs.*ix;
py = meany + hcs.*iy;
qy = meany - hcs.*iy;
pix = coord2ind(GRID,px,py);
qix = coord2ind(GRID,qx,qy);

% get values
p = nan(size(pix));
q = p;
nx = ~isnan(pix) & ~isnan(qix);
p(nx) = GRID.Z(pix(nx));
q(nx) = GRID.Z(qix(nx));

if nargout==1
    if nargin==3
        fct = validatestring(varargin{1},...
            {'min','max','mean','diff','normdiff'});
    else
        fct = 'mean';
    end
    switch fct
        case 'mean'
            val = nanmean([p,q],2);
        case 'max'
            val = nanmax([p,q],[],2);
        case 'min'
            val = nanmin([p,q],[],2);
        case 'diff'
            val = abs(diff([p,q],1,2));
        case 'normdiff'
            val = abs(diff([p,q],1,2))./nansum([p,q],2);
    end
    varargout{1} = val;
else
    varargout{1} = p;
    varargout{2} = q;
end
if nargout>2
    varargout{3} = pix;
    varargout{4} = qix;
end


end

