function [Z,R] = GRIDobj2im(DEM)

% GRIDOBJ2IM Create image from GRIDobj
%
% Syntax
%
%      [Z,R] = GRIDobj2im(DEM)
%
% Description
%
%      GRIDobj2im creates an image Z and an imref2d object R from a GRIDobj
%      DEM. The matrix Z is flipped along the first dimension so that
%      Z(1,1) refers to the southwest corner. This is necessary to have a
%      valid imref2d object.
%      
% Input arguments
%   
%      DEM     GRIDobj
%
% Output arguments
%
%      Z       grayvalue matrix
%      R       imref2d reference to world coordinates
%
% Example 1
%
%      DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%      [Z,R] = GRIDobj2im(DEM);
%      imshow(Z,R,[]);
%      axis xy
%      colormap(parula)
%
% Example 2
%
%      DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%      [~,R] = GRIDobj2im(DEM);
%      RGB = imageschs(DEM,[],'colormap','landcolor');
%      % Following line is required because columns must start south
%      % but GRIDobj standard is columns-start-north
%      RGB = flipud(RGB); 
%      imshow(RGB,R,[]);
%      axis xy
%
% See also: FLOWobj, STREAMobj, GRIDobj/info
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 17. August, 2017

[xWorldLimits,yWorldLimits] = getcoordinates(DEM);
xWorldLimits = xWorldLimits([1 end]);
yWorldLimits = yWorldLimits([1 end]);

xWorldLimits = xWorldLimits(:)';
yWorldLimits = yWorldLimits(:)';
if yWorldLimits(1) > yWorldLimits(2)
    yWorldLimits = fliplr(yWorldLimits);
    Z = flipud(DEM.Z);
else
    Z = DEM.Z;
end

R = imref2d(DEM.size,xWorldLimits,yWorldLimits);