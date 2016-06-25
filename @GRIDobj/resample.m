function DEMr = resample(DEM,target,method)

% resample grid to alter spatial resolution
%
% Syntax
%
%     DEMr = resample(DEM,cellsize)
%     DEMr = resample(DEM,GRID)
%     DEMr = resample(...,'method')
%
% Description
%
%     resample changes the cellsize of a grid. The function uses the Matlab
%     function imtransform. If an instance of GRIDobj is supplied as
%     second argument, resample interpolates values in DEM to match the
%     spatial reference of GRID.
%
% Input arguments
%
%     DEM       grid object (GRIDobj)
%     cellsize  cellsize of resampled grid
%     GRID      other grid object
%     method    'bicubic', 'bilinear', or 'nearest' 
%
% Output arguments
%
%     DEMr    grid object (GRIDobj)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     DEMr = resample(DEM,5);
%     imagesc(DEMr)
%
% See also: griddedInterpolant
%
%        
% Author:  Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 8. August, 2015 

% check input arguments
narginchk(2,3)
validateattributes(target,{'double' 'GRIDobj'},{'scalar'})
if nargin == 2;
    method = 'bilinear';
else
    method = validatestring(method,{'bicubic', 'bilinear', 'nearest' });
end

% check underlying class
if islogical(DEM.Z)
    method = 'nearest';
end

% get coordinate vectors
[u,v] = getcoordinates(DEM);

% tform
T = maketform('affine',[1 0 0; 0 1 0; 0 0 1]);

% Fillvalues
if isinteger(DEM.Z)
    fillval = 0;
elseif islogical(DEM.Z)
    fillval = 0;
else
    fillval = nan;
end


if isa(target,'GRIDobj')
    
    DEMr    = target;
    [xn,yn] = getcoordinates(DEMr);
    
    DEMr.Z = imtransform(DEM.Z,T,method,...
        'Udata',[u(1) u(end)],'Vdata',[v(1) v(end)],...
        'Xdata',[xn(1) xn(end)],'Ydata',[yn(1) yn(end)],...
        'Size',DEMr.size,...
        'FillValues',fillval);
    DEMr.name = [DEM.name ' (resampled)'];
        
else
    csnew   = target;
    DEMr    = GRIDobj([]);
    [DEMr.Z,xn,yn] = imtransform(DEM.Z,T,method,...
        'Udata',[u(1) u(end)],'Vdata',[v(1) v(end)],...
        'Xdata',[u(1) u(end)],'Ydata',[v(1) v(end)],...
        'XYscale',[csnew csnew],...
        'FillValues',fillval);


    % new referencing matrix
    DEMr.refmat = [0 -csnew;...
                   csnew 0; ...
                   xn(1)-csnew ...
                   yn(1)+csnew];
    % size of the resampled grid           
    DEMr.size    = size(DEMr.Z);
    DEMr.cellsize = csnew;
    DEMr.georef = DEM.georef;
end

DEMr.name    = [DEM.name ' (resampled)'];

if ~isempty(DEMr.georef) && ~isa(target,'GRIDobj');
    DEMr.georef.RefMatrix = DEMr.refmat;
    DEMr.georef.Height = DEMr.size(1);
    DEMr.georef.Width  = DEMr.size(2);
    
    DEMr.georef.SpatialRef = refmatToMapRasterReference(DEMr.refmat, DEMr.size);
    
end








