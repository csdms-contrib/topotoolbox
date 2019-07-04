function OUT = tidy(SW,varargin)

%TIDY remove overlapping points from SWATHobj
%
% Syntax
%
%     OUT = tidy(SW)
%     OUT = tidy(SW,I)
%     OUT = tidy(SW,I,'both')
%
% Description
%
%     TIDY(SW) removes overlapping points in a SWATHobj by setting the
%     corresponding Z values to NaN. 
%
%
% Input arguments
%
%     SW     instance of SWATHobj
%     I      instance of GRIDobj
%
% Output arguments
%
%     OUT    instance of SWATHobj
%
% Examples
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     FA  = flowacc(FD);
%     S = STREAMobj(FD,FA>1000);
%     S = trunk(klargestconncomps(S,1));
%     [x,y] = STREAMobj2XY(S);
%     ix = ~isnan(x);
%     SW = SWATHobj(DEM,x(ix),y(ix),'smooth',1e3,'width',3e3);
%     SW = tidy(SW);
%     figure,plot(SW,'points',true),axis equal
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: March, 2017


% check inputs
narginchk(1,3)

if ~isa(SW,'SWATHobj')
    error('First input needs to be of class SWATHobj')
end

if nargin>=2
    I = varargin{1};
    if ~isa(FD,'GRIDobj')
        error('Second input needs to be of class GRIDobj')
    end
end
if nargin==3
    if ~strcmp('both',varargin{2})
        error('Unknown third input to function call')
    end
end


% Output variable
OUT = SW;

if nargin==1 || nargin==3
    
    % Use the profile points that define the swath
    xs = SW.xy(:,1);
    ys = SW.xy(:,2);
    
    % Add points at great distance to prevent polygons with 'inf'
    if (1)
        theta = 0:5:359;
        rho = 1e10.*ones(size(theta)); % great distance here = 1e10 map units
        [dx,dy] = pol2cart(theta,rho);
        edgex = mean(xs)+dx;
        edgey = mean(ys)+dy;
        xs = [xs;edgex'];
        ys = [ys;edgey'];
    end
    
    % Triangulate and create voronoi polygons
    dt = DelaunayTri(xs(:),ys(:));
    [V,C] = voronoiDiagram(dt);
    
    % Loop over voronoi polygons - one polygon per point along profile
    IX = zeros(size(SW.Z));
    for i = 1 : length(SW.Z(1,:))
        xp = SW.X(:,i);
        yp = SW.Y(:,i);
        ix = C{i};
        xv = V(ix,1);
        yv = V(ix,2);
        % determine which points are inside the polygon
        IN = inpolygon(xp,yp,xv,yv);
        IX(:,i) = IN;
    end
    % blank out points outside the polygons
    OUT.Z(not(IX)) = nan;
    
end



if exist('I','var')
    
    
    
    
    
end



end



