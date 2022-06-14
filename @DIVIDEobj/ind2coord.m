function [x,y] = ind2coord(D,ix)
%IND2COORD convert linear index to x and y coordinates
%
% Syntax
%
%     [x,y] = ind2coord(D,ix)
%
% Description
%
%     ind2coord converts a linear index into an instance of DIVIDEobj to x
%     and y coordinates. Note that the virtual grid associated
%     with a DIVIDEobj is shifted by half a cell size in x and y and has an
%     extra row and collumn. The linear indices point to the coreners of
%     the pixels in the GRIDobj the DIVIDEobj is based on.
%
% Input arguments
%
%     D       instance of DIVIDEobj
%     ix      linear index
%
% Output arguments
%
%     x,y     x- and y-coordinates  
%
% See also: GRIDobj/ind2coord, GRIDobj/coord2ind, GRIDobj/getcoordinates
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
%         code simply copied and help adjusted by Dirk Scherler
% Date: 30. January, 2013 and February 2019


[r,c] = ind2sub(D.size,ix(:));
xy    = double([r c ones(numel(ix),1)])*D.refmat;
x = xy(:,1);
y = xy(:,2);