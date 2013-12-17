function DEMr = resample(DEM,target,method)

% resample grid to alter spatial resolution
%
% Syntax
%
%     DEMr = resample(DEM,scale)
%     DEMr = resample(DEM,GRID)
%     DEMr = resample(...,'method')
%
% Description
%
%     resample changes the cellsize of a grid. The function uses the Matlab
%     class griddedInterpolant. If an instance of GRIDobj is supplied as
%     second argument, resample interpolates values in DEM to match the
%     spatial reference of GRID.
%
% Input arguments
%
%     DEM     grid object (GRIDobj)
%     scale   resampling scale. A scale of 2 will change the cellsize x to
%             x/2. Scales larger/less than 1 will increase/decrease 
%             resolution.
%     GRID    other grid object
%     method  'linear' (default). See griddedInterpolant for more options
%
% Output arguments
%
%     DEMr    grid object (GRIDobj)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     DEMr = resample(DEM,.2);
%     imagesc(DEMr)
%
% See also: griddedInterpolant
%

% check input arguments
narginchk(2,3)
validateattributes(target,{'double' 'GRIDobj'},{'scalar'})
if nargin == 2;
    method = 'linear';
else
    method = validatestring(method,{'linear','nearest','spline','pchip','cubic'});
end

% get coordinate vectors
[x,y] = refmat2XY(DEM.refmat,DEM.size);

if isa(target,'GRIDobj')
    DEMr    = target;
    [xn,yn] = refmat2XY(DEMr.refmat,DEMr.size);
    if yn(2)<yn(1)
        yn = yn(end:-1:1);
    end
else
    scale   = target;
    DEMr    = GRIDobj([]);
    DEMr.georef = DEM.georef;
    % cellsize of resampled grid
    cellsizenew = DEM.cellsize/scale;
    % coordinate vectors
    xn    = min(x):cellsizenew:max(x);
    yn    = min(y):cellsizenew:max(y);

    % new referencing matrix
    DEMr.refmat = [0 -cellsizenew;...
                   cellsizenew 0; ...
                   xn(1)-cellsizenew ...
                   yn(end)+cellsizenew];
    % size of the resampled grid           
    DEMr.size    = [numel(yn) numel(xn)];
    DEMr.cellsize = cellsizenew;
end
DEMr.name    = [DEM.name ' (resampled)'];

% griddedInterpolant requires coordinate vectors monotically increasing
if y(2)<y(1)
    y = y(end:-1:1);
    Z = flipud(DEM.Z);
else
    Z = DEM.Z;
end
% gridded interpolant
cla   = class(Z);
F     = griddedInterpolant({x(:)' y(:)},double(Z'),method);

% interpolate
DEMr.Z = flipud(F({xn(:)' yn(:)})');

if islogical(Z)
    DEMr.Z(isnan(DEMr.Z)) = 0;
end
    
DEMr.Z = cast(DEMr.Z,cla);


if ~isempty(DEMr.georef)
    DEMr.georef.RefMatrix = DEMr.refmat;
    DEMr.georef.Height = DEMr.size(1);
    DEMr.georef.Width  = DEMr.size(2);
    
    DEMr.georef.SpatialRef = refmatToMapRasterReference(DEMr.refmat, DEMr.size);
    
end








