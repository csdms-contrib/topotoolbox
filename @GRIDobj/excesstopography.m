function DEM = excesstopography(DEM,varargin)

% difference between actual elevations and elevations with threshold slope
% 
% Syntax
%
%    EXT = excesstopography(DEM)
%    EXT = excesstopography(DEM,pn,pv,...)
%
% Description
%
%    excesstopography uses grayscale morphological erosion to reconstruct
%    an idealised surface that features only slopes equal or less than the 
%    threshold slope. Subsequently, the function calculates the difference
%    between the DEM and idealized topography which is referred to excess
%    topography.
%
%    Note that this function takes a while to evaluate, especially if the
%    DEM and kernelsize are large. Be sure that the kernel is large enough.
%
% Input arguments
%
%    DEM    Digital elevation model (class: GRIDobj)
%    
% Parameter name/value pairs
%
%    'maxgradient'   maximum gradient in degrees of slopes (default: 30°)
%    'kernelsize'    side length in pixels used by the kernel. Must be 
%                    integer and odd (default: 101)
%    'output'        'difference' (default) or 'elevation'. Latter returns
%                    the eroded DEM without calculating the difference.
%
% Output arguments
%
%    EXT    excess topography (class: GRIDobj)
%
% Example
%
%    DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%    EXT = excesstopography(DEM,'maxgradient',30,'kernelsize',31);
%    imageschs(DEM,EXT)
%  
% Reference
%
%    Blöthe, J.H., Korup, O., Schwanghart, W. (2015): Large landslides lie low: 
%    Excess topography in the Himalaya-Karakorum ranges. Geology, 43, 523-526. 
%    [DOI: 10.1130/G36527.1] 
%
% See also: GRIDobj/localtopography 
%     
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 29. April, 2015



% Parse Inputs
p = inputParser;         
p.FunctionName = 'excesstopography';
addRequired(p,'DEM', @(x) isa(x,'GRIDobj'));
addParamValue(p,'maxgradient',30,@(x) isscalar(x) && (x >= 0));
addParamValue(p,'kernelsize',101, @(x) isscalar(x) && (rem(x,2)==1) && (x>=3));
addParamValue(p,'output','difference', @(x) ischar(validatestring(x,{'difference','elevation'})));

parse(p,DEM,varargin{:});

maxgradient = p.Results.maxgradient;
kernelsize  = p.Results.kernelsize;
output      = validatestring(p.Results.output,{'difference','elevation'});

cl          = class(DEM.Z);

% maxgradient is provided in degrees. Convert to radians
maxgradient = tand(maxgradient);

% prepare for ordfilt
% handle nans and edge effects
INAN   = isnan(DEM.Z);
m      = max(DEM);
DEM.Z(INAN) = inf;

% prepare kernel and offset
% rectangular kernel
domain = ones(kernelsize);
% offset 
offset = false(size(domain));
center = ceil(kernelsize/2);
offset(center,center) = true;
offset = bwdist(offset,'e');
% offset = max(offset(:))-offset;
offset = offset * DEM.cellsize * maxgradient;

% do the calculation
switch output
    case 'difference'
        DEM   = DEM - (ordfilt2(double(DEM.Z)-m,1,domain,double(offset),'zeros') + m);
        DEM.Z = cast(DEM.Z,cl);
    case 'elevation'
        DEM.Z = ordfilt2(double(DEM.Z)-m,1,domain,double(offset),'zeros') + m;
        DEM.Z = cast(DEM.Z,cl);
end
% replace infs with nan again
DEM.Z(INAN) = nan;
