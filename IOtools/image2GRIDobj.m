function DEM = image2GRIDobj(I,R)

%IMAGE2GRIDOBJ convert image to GRIDobj
%
% Syntax
%
%     A = image2GRIDobj(filename,cellsize)
%     A = image2GRIDobj(I,cellsize)
%     A = image2GRIDobj(...,B)
%
% Description
%
%     image2GRIDobj reads a single band image and converts it to a GRIDobj
%     with a spatial reference defined by the cellsize.
%
% Input arguments
%
%     I          grayscale image
%     cellsize   cellsize
%     B          GRIDobj
%
% Output arguments
%
%     A        GRIDobj
%
%
% See also: GRIDobj, imread
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 15. March, 2018


if ischar(I) || isstring(I)
    name = I;
    I = imread(I);
    siz = size(I);
end

% create

if isa(R,'map.rasterref.MapCellsReference')
    x = R.
    
