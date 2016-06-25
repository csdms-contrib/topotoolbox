function [varargout] = SWATHobj2GRIDobj(SW,DEM,varargin)
% SWATHOBJ2GRIDOBJ creates a GRIDobj with Swath-specific information
%
% Syntax
%
%     OUT = SWATHobj2GRIDobj(SW,DEM)
%     OUT = SWATHobj2GRIDobj(SW,DEM,type)
%     [OUT1,OUT2,...] = SWATHobj2GRIDobj(SW,DEM,type1,type2,...)
%
% Description
%
%     SWATHobj2GRIDobj(SW,DEM)
%
%
%     Make sure that the SW is fully contained withtin the DEM
% 
%
% Input arguments
%
%     SW     instance of SWATHobj
%     DEM    instance of GRIDobj
%     type   string that specifies the type of output. Valid choices are
%            'distx', 'disty', and 'ix'. 'distx' and 'disty' give the 
%            distance along and and across the swath, respectively. 'ix' 
%            creates a binary image indicating pixels inside the SWATHobj.
%
% Output arguments
%
%     OUT    instance of GRIDobj
%
% Examples
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     A  = flowacc(FD);
%     S = STREAMobj(FD,A>100);
%     S = trunk(klargestconncomps(S,1));
%     [x,y] = STREAMobj2XY(S);
%     ix = ~isnan(x);
%     xval = flipud(x(ix));
%     yval = flipud(y(ix));
%     SW = SWATHobj(DEM,xval,yval,'smooth',200);
%     %   SW = tidy(SW);
%     DX = SWATHobj2GRIDobj(SW,DEM,'distx');
%     figure,imagesc(DX), hold on
%     plot(SW,'points',false,'labeldist',1e3), hold off
%     colorbar
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: May, 2015


narginchk(2,5)

if length(varargin)~=nargout
    error('Number of inputs mismatches number of outputs.')
end

if nargin == 2;
    type{1} = 'ix';
else
    for i = 1 : length(varargin)
        validatestring(varargin{i},{'distx','disty','ix','z'});
    end
    type = varargin;
end


[x,y] = getcoordinates(DEM);
[X,Y] = meshgrid(x,y);

% Delineate SWATHobj in GRID domain
% (might be better to use inpolygon(X,Y,xv,yv)?)
IM = true(size(SW.Z));
IM(isnan(SW.Z)) = false;
B = bwboundaries(IM,4);
A = false(DEM.size);
for k = 1 : length(B)
    ix_bdy = sub2ind(size(IM),B{k}(:,1),B{k}(:,2));
    x_bdy = SW.X(ix_bdy);
    y_bdy = SW.Y(ix_bdy);
    if length(x_bdy)>2
        [newx,newy] = eqdistline(x_bdy,y_bdy,DEM.cellsize/3);
    else
        newx = x_bdy; newy = y_bdy;
    end
    ix_edge = coord2ind(DEM,newx,newy);
    ix_edge = unique(ix_edge);
    ix_edge = ix_edge(~isnan(ix_edge));
    A(ix_edge) = true;
end
A = imfill(A,4,'holes');

% Swath locations
ix_SW = ~isnan(SW.Z);
x = SW.X(ix_SW);
y = SW.Y(ix_SW);

% Grid locations
ix_grid = find(A==1);
xq = X(ix_grid);
yq = Y(ix_grid);
[xq_edges,yq_edges] = ind2coord(DEM,ix_edge);

varargout = cell(size(type));
for i = 1 : length(type)
    if strcmp(type{i},'ix')
        OUT = DEM;
        OUT.Z(:) = false;
        OUT.Z(ix_grid) = true;
        OUT.Z(ix_edge) = true;
        varargout{i} = OUT;
    else
        switch type{i}
            case 'distx' % longitudinal distance
                D = repmat(SW.distx',length(SW.disty),1);
            case 'disty' % transverse distance
                D = repmat(SW.disty,1,length(SW.distx));
            case 'z' % z-values of the SWATHobj
                D = SW.Z;
        end
        OUT = DEM;
        OUT.Z(:) = nan;
        % interpolate between given points
        v = D(ix_SW);
        F = TriScatteredInterp(x,y,v);
        vq = F(xq,yq);
        OUT.Z(ix_grid) = vq;
        % the edges might be a bit problematic
        F = TriScatteredInterp(x,y,v,'nearest');
        vq = F(xq_edges,yq_edges);
        OUT.Z(ix_edge) = vq;
        varargout{i} = OUT;
    end
end

