function [U,V] = M2UV(M,X,Y,px,W)

% calculate horizontal, directional components from flow direction matrix
%
% Syntax
%
%     [U,V] = M2UV(M,X,Y,lag,W)
%
% Description
%
%     M2UV calculates the horizontal flow vectors of a digital 
%     elevation model based on the flow direction matrix. The vector
%     quantities in the X and Y directions are calculated as weighted
%     average from the multiple flow direction matrix.
%
% Input
%
%     M         flow direction matrix as returned by flowdir
%     X,Y       coordinate matrices
%     lag       spatial lag (smoothes the vector field (default = 1), but
%               removes values on edges)
%     W         velocity matrix same size as X. Scales the directional
%               components U and V (U.*W, V.*W).
% 
% Output
%
%     U,V       directional components
%                  
% Example
% 
%     load exampleDEM
%     M = flowdir(X,Y,fillsinks(dem));
%     A = flowacc(M,size(dem));
%     [U,V] = M2UV(M,X,Y,2,log(A));
%     theta = cart2pol(U,V);
%     imageschs(X,Y,dem,theta);
%     hold on
%     quiver(X,Y,U,V,'k')
%
%
% See also: QUIVER, CART2POL
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 22. July, 2010



siz = size(X);


if nargin >= 4;
    M   = M^px;
    nrc = numel(X);
    M   = spdiags(spfun(@(x) 1./x,sum(M,2)),0,nrc,nrc) * M;
end

% identify sink pixels and pixels on edge
I = reshape(sum(M,2) == 0,siz);
V = reshape(M*Y(:),siz)-Y;
U = reshape(M*X(:),siz)-X;

U(I) = nan;
V(I) = nan;

if nargin == 5;
    U = U.*W;
    V = V.*W;
end