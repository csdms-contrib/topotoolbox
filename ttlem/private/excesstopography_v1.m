function DEM = excesstopography_v1(DEM,varargin)

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
%    DEM and kernelsize are large. Be sure that the kernel is large enough
%    or set the option 'iterate' to true.
%
% Input arguments
%
%    DEM    Digital elevation model (class: GRIDobj)
%    
% Parameter name/value pairs
%
%    'maxgradient'   maximum gradient in degrees of slopes (default: 30°)
%    'unit'          deg (degrees = default) or tan (tangens). Applies to
%                    'maxgradient' and 'tol'.
%    'kernelsize'    side length in pixels used by the kernel. Must be 
%                    integer and odd (default: 7)
%    'output'        'difference' (default) or 'elevation'. Latter returns
%                    the eroded DEM without calculating the difference.
%    'iterate'       true (default) or false. A small kernel may not
%                    detect all excesstopography. Setting iterate to true
%                    repeatedly erodes the topography until maxgradient is
%                    reached
%    'tol'           only applicable if 'iterate' is set to true. Iteration
%                    terminates if max(gradient8(DEM))-maxgradient < tol.
%                    The default is 6e-4 ° (degrees). Note that setting tol
%                    to a very small value might cause that the algorithm
%                    uses all maxiter iterations.
%    'maxiter'       only applicable if 'iterate' is set to true. Iteration
%                    terminates if maxiter are reached. The default is 100.
%                    
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
%
% See also: GRIDobj/localtopography 
%     
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 23. January, 2014



% Parse Inputs
p = inputParser;         
p.FunctionName = 'excesstopography';
addRequired(p,'DEM', @(x) isa(x,'GRIDobj'));
addParamValue(p,'maxgradient',30,@(x) isscalar(x) && (x >= 0));
addParamValue(p,'kernelsize',101, @(x) isscalar(x) && (rem(x,2)==1) && (x>=3));
addParamValue(p,'output','difference', @(x) ischar(validatestring(x,{'difference','elevation'})));
addParamValue(p,'iterate',true);
addParamValue(p,'maxiter',100);
addParamValue(p,'tol',6e-4);
addParamValue(p,'unit','deg', @(x) ischar(validatestring(x,{'deg','tan'})));

parse(p,DEM,varargin{:});

maxgradient = p.Results.maxgradient;
kernelsize  = p.Results.kernelsize;
output      = validatestring(p.Results.output,{'difference','elevation'});
iterate     = p.Results.iterate;
tol         = p.Results.tol;
unit        = validatestring(p.Results.unit,{'deg','tan'});

cl          = class(DEM.Z);

% maxgradient is provided in degrees. Convert to radians
switch unit
    case 'deg'
        maxgradient = tand(maxgradient);
        tol = tand(tol);
end

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
offset_ori = bwdist(offset,'e');
% offset = max(offset(:))-offset

% do the calculation
DEMcopy = DEM;
%%
G  = gradient8(DEMcopy);
mG=max(G);
G_s=GRIDobj(DEM);
[dx,dy]=gradient(DEMcopy.Z,DEM.cellsize);
G_s.Z=sqrt(dx.^2+dy.^2);
mG_s=max(G_s);
mG_s=max(mG_s,mG);
maxgradientIni=maxgradient;
iter_m =0;

while round(mG_s*100)>round(maxgradientIni*97.5) && (iter_m<= p.Results.maxiter);
    iter_m = iter_m+1;      
    if mG > maxgradient;
        offset = offset_ori * DEM.cellsize * maxgradient;
        if iterate
            iter = 0;
            while ((mG-maxgradient) > tol) && (iter<= p.Results.maxiter);
                iter = iter+1;                
                DEMcopy.Z = ordfilt2(double(DEMcopy.Z)-m,1,domain,double(offset),'zeros') + m;
                G  = gradient8(DEMcopy);
                mG = max(G);
            end
        else
            DEMcopy.Z = ordfilt2(double(DEMcopy.Z)-m,1,domain,double(offset),'zeros') + m;
        end
    end
    [dx,dy]=gradient(DEMcopy.Z,DEM.cellsize);
    G_s.Z=sqrt(dx.^2+dy.^2);
    mG_s=max(G_s);
    maxgradient=0.975*maxgradient;
end

%%
switch output
    case 'difference'
        DEM   = DEM-DEMcopy;
        DEM.Z = cast(DEM.Z,cl);
    case 'elevation'
        DEM.Z = cast(DEMcopy.Z,cl);
end
% replace infs with nan again
DEM.Z(INAN) = nan;
DEM.name = 'excess topography';