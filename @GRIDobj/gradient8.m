function G = gradient8(DEM,unit)

% 8-connected neighborhood gradient and aspect of a digital elevation model
%
% Syntax
%
%     G = gradient8(DEM,unit)
%
% Description
%
%     gradient8 returns the numerical steepest downward gradient and aspect 
%     of a digital elevation model using an 8-connected neighborhood. 
%
% Input
%
%     DEM       digital elevation model (class: GRIDobj)
%     unit      'tan' --> tangent (default)
%               'rad' --> radian
%               'deg' --> degree
%               'sin' --> sine
%               'per' --> percent
% 
% Output
%
%     G         gradient (class: GRIDobj)
%                  
% Example
% 
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     G = gradient8(DEM,'degree');
%     subplot(2,1,1)
%     imagesc(DEM)
%     subplot(2,1,2)
%     imagesc(G)
%
%
% See also: GRIDobj, CURVATURE, ASPECT
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 6. November, 2015

if nargin == 1;
    unit = 'tangent';
else
    unit = validatestring(unit,{'tangent' 'degree' 'radian' 'percent' 'sine'},'gradient8','unit',2);
end

% create a copy of the DEM instance
G = DEM;
c = class(DEM.Z);
switch c
    case 'double'
        G.Z = double.empty(0,0);
    otherwise
        G.Z = single.empty(0,0);
        c   = 'single';
end

% I found Large matrix support using blockproc inefficient for gradient8.
% Matrix dimensions have thus been increased to an out-of-range value to
% avoid calling blockproc.
% Large matrix support. Break calculations in chunks using blockproc.

if numel(DEM.Z)>(20001*20001);
    blksiz = bestblk(size(DEM.Z),5000);
    c   = class(DEM.Z);
    
    switch c
        case {'double', 'single'}
            padval = inf;
        case 'logical'
            padval = true;
        otherwise
            padval = intmax(c);
    end
    cs  = G.cellsize;
    fun = @(x) steepestgradient(x,cs,c);
    G.Z = blockproc(DEM.Z,blksiz,fun,...
           'BorderSize',[1 1],...
           'Padmethod',padval,...
           'UseParallel',true);
else
    G.Z = steepestgradient(DEM.Z,G.cellsize,c);
end

G.name = 'gradient';
G.zunit = unit;

switch unit
    case 'tangent'
        % do nothing
    case 'degree'
        G.Z = atand(G.Z);
    case 'radian'
        G.Z = atan(G.Z);
    case 'sine'
        G.Z = sin(atan(G.Z));
    case 'percent'
        G.Z = G.Z*100;
end
end




function G = steepestgradient(z,cellsize,c)

if isstruct(z);
    z = z.data;
end
    

% check for nans;
I = isnan(z);
flagnan = any(I(:));
if flagnan
    z(I) = inf;
end

NEIGH = false(3);
% calculate along orthogonal neighbors
NEIGH(2,:) = true;
NEIGH(:,2) = true;

switch c
    case 'double'
        G = (z-imerode(z,NEIGH))/cellsize;
    case 'single'
        G = single(z-imerode(z,NEIGH))/cellsize;
end

% calculate along diagonal neighbors
NEIGH(:,:) = false;
NEIGH([1 5 9 3 7]) = true;

switch c
    case 'double'
        G = max(G,(z-imerode(z,NEIGH))/norm([cellsize cellsize]));
    case 'single'
        G = max(G,single(z-imerode(z,NEIGH))/single(norm([cellsize cellsize])));
end

if flagnan
    G(I) = nan;
end
end

