function OUT2 = hillshade(DEM,varargin)

% create hillshading from a digital elevation model (GRIDobj)
%
% Syntax
%    
%     H = hillshade(DEM)
%     H = hillshade(DEM,'pn','pv',...)
%
% Description
%
%     Hillshading is a very powerful tool for relief depiction.
%     hillshade calculates a shaded relief for a digital elevation model 
%     based on the angle between the surface and the incoming light beams.
%     If no output arguments are defined, the hillshade matrix will be
%     plotted with a gray colormap. The hillshading algorithm follows the
%     logarithmic approach to shaded relief representation of Katzil and
%     Doytsher (2003).
%
% Input
%
%     DEM       Digital elevation model (class: GRIDobj)
%
% Parameter name/value pairs
%
%     azimuth        azimuth angle, (default=315)
%     altitude       altitude angle, (default=60)
%     exaggerate     elevation exaggeration (default=1). Increase to
%                    pronounce elevation differences in flat terrain
%
% Output
%
%     H         shaded relief (ranges between 0 and 1)
%
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     hillshade(DEM)
% 
% References
%
%     Katzil, Y., Doytsher, Y. (2003): A logarithmic and sub-pixel approach
%     to shaded relief representation. Computers & Geosciences, 29,
%     1137-1142.
%
% See also: SURFNORM, IMAGESCHS
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 21. December, 2012



% Parse inputs
p = inputParser;
p.StructExpand  = true;
p.KeepUnmatched = false;
p.FunctionName = 'hillshade'; 
addRequired(p,'DEM',@(x) isequal(class(x),'GRIDobj'));
addParamValue(p,'azimuth',315,@(x) isscalar(x) && x>= 0 && x<=360);
addParamValue(p,'altitude',60,@(x) isscalar(x) && x>= 0 && x<=90);
addParamValue(p,'exaggerate',1,@(x) isscalar(x) && x>0);
parse(p,DEM,varargin{:});

OUT     = DEM;
OUT.Z   = [];

% Large matrix support. Break calculations in chunks using blockproc
if numel(DEM.Z)>(1001*1001);
    blksiz = bestblk(size(DEM.Z),1000);    
    padval = 'symmetric';
    OUT.Z = blockproc(DEM.Z,blksiz,@(x) hsfun(x),...
                'BorderSize',[1 1],...
                'padmethod',padval,...
                'UseParallel',true);
else
    OUT.Z = hsfun(DEM.Z);
end

OUT.name = 'hillshade';
OUT.zunit = '';

if nargout == 0;
    OUT.Z = uint8(OUT.Z*255);
    imagesc(OUT);
    colormap(gray)
else
    OUT2 = OUT;
end


%% Subfunction
function H = hsfun(Z)

if isstruct(Z)
    Z = Z.data;    
end

% correct azimuth so that angles go clockwise from top
azid = p.Results.azimuth-90;

% use radians
altsource = p.Results.altitude/180*pi;
azisource = azid/180*pi;

% calculate solar vector
[sx,sy,sz] = sph2cart(azisource,altsource,1);

% calculate surface normals
[Nx,Ny,Nz] = surfnorm(Z/DEM.cellsize*p.Results.exaggerate);

% calculate cos(angle)
H = [Nx(:) Ny(:) Nz(:)]*[sx;sy;sz];

% reshape
H = reshape(H,size(Nx)); 

% % usual GIS approach
% H = acos(H);
% % force H to range between 0 and 1
% H = H-min(H(:));
% H = H/max(H(:));


end
end