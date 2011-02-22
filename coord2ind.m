function [IX,ixcoord,res] = coord2ind(X,Y,x,y)

% convert xy coordinates to linear index
%
% Syntax
%
%     IX = coord2ind(X,Y,x,y)
%
% Description
%
%     coord2ind converts xy coordinates to a linear index. 
%
% Input
%
%     X,Y       coordinate matrices. X and Y must be plaid and
%               monotonically increasing as produced by meshgrid
%     x,y       point coordinates
%
% Output
%
%     IX        vector with linear indices of the point coordinates in
%               X and Y
%     ixcoord   linear index such that 
%               IX = coord2ind(X,Y,x(ixcoord),y(ixcoord))
%     res       distance between original x y data and clipped data 
%
% Example
%
%     [X,Y] = meshgrid(0:0.01:1,0:0.01:1);
%     x = rand(100,1);
%     y = rand(100,1);
%     Z = false(size(X));
%     [IX,ix,res] = coord2ind(X,Y,x,y);
%     Z(IX) = true;
%     imagesc(X(1,:),Y(:,2),Z); axis image; axis xy
%     hold on
%     scatter(x,y,'y') 
%
% See also: SUB2IND, IND2SUB, MESHGRID
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 15. March, 2009



% force column vectors
x = x(:);
y = y(:);

siz = size(X);

dx = X(1,2)-X(1);
dy = Y(2)-Y(1);

IX1 = (x-X(1))./dx + 1;
IX2 = (y-Y(1))./dy + 1;

IX1 = round(IX1);
IX2 = round(IX2);

I = IX1>siz(2) | IX1<1 | IX2>siz(1) | IX2<1;

x(I) = [];
y(I) = [];
ixcoord = find(~I);

if ~isempty(x)
    IX = sub2ind(siz,IX2,IX1);
    if nargout == 3
        res = hypot(X(IX)-x,Y(IX)-y);
    end
else
    IX = [];
    ixcoord = [];
    res = [];
end




