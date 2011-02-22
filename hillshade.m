function h = hillshade(X,Y,dem,azid,altd,exagg)

% create hillshading from a digital elevation model
%
% Syntax
%    
%     H = hillshade(X,Y,dem,azid,altd,exagg)
%     hillshade(X,Y,dem,azid,altd,exagg)
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
%     X,Y       coordinate matrices as returned by meshgrid
%     dem       digital elevation model
%     azid      azimuth angle in degrees (default: 315)
%     altd      altitude angle in degrees (default: 60)
%     exagg     elevation exaggeration (positive scalar, default: 1)
%
% Output
%
%     H         shaded relief (ranges between 0 and 1)
%
%
% Example
%
%     load exampleDEM
%     H = hillshade(X,Y,dem,315,40);
%     imagesc(X(1,:),Y(:,2),H); axis image; axis xy
%     colormap('gray')
% 
% References
%
%     Katzil, Y., Doytsher, Y. (2003): A logarithmic and sub-pixel approach
%     to shaded relief representation. Computers & Geosciences, 29,
%     1137-1142.
%
% See also: SURFNORM, IMAGESCHS
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 18. May, 2009

error(nargchk(3, 6, nargin))

defazid = 315;
defaltd = 60;
defexagg= 1;

if nargin == 3
    azid   = defazid;
    altd   = defaltd;
    exagg  = defexagg;
elseif nargin == 4;
    altd   = defaltd;
    exagg  = defexagg;
elseif nargin == 5;
    exagg  = defexagg;
end

if isempty(azid); azid  = defazid; end
if isempty(altd); altd  = defaltd; end
if isempty(exagg);exagg = defexagg;end

if Y(1)>Y(2);
    flipud(Y);
    flipud(dem);
    flipd = true;
else
    flipd = false;
end



if ~isscalar(azid) || ~isscalar(altd)
    error('TopoToolbox:incorrectinput',...
          'azid and thetad must be scalars')
end

siz  = size(dem);

% correct azimuth so that angles go clockwise from top
azid = -azid+90;

% use radians
altsource = altd/180*pi;
azisource = azid/180*pi;

% calculate surface normals
[Nx,Ny,Nz] = surfnorm(X,Y,dem*exagg);

% calculate solar vector
[sx,sy,sz] = sph2cart(azisource,altsource,1);

% calculate cos(angle)
H = [Nx(:) Ny(:) Nz(:)]*[sx;sy;sz];

% reshape
H = reshape(H,siz); 

% % usual GIS approach
% H = acos(H);
% % force H to range between 0 and 1
% H = H-min(H(:));
% H = H/max(H(:));

% logarithmic approach
H = 1/log(2) * log(1+H.^4);
 
if ~flipd
    H = H*(-1)+1;
end

% check output arguments
if nargout == 0
    imagesc(X(1,:),Y(:,2),H); 
    axis image; 
    axis xy;
    colormap(flipud(gray));
else
    h=H;
end


