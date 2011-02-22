function I = dependencemap(M,I)

% drainage area for specific locations in a digital elevation model
%
% Syntax
%
%     D = dependencemap(M,L)
%
% Description
%
%     dependence map returns a logical matrix masking the upslope part
%     of the digital elevation model that drains into specified cells.
%
% Input
%
%     M         multiple or single flow direction matrix
%     L         logical matrix same size as the digital elevation from
%               which the flow direction matrix was derived from. true
%               values in L indicate cells whose upslope area is to be
%               calculated
%
% Output
%
%     D         logical dependence matrix
%  
% Example
%
%     load exampleDEM
%     M = flowdir(X,Y,dem,'exponent',5);
%     M = multiremfracs(M,0.05);
%     L = false(size(dem));
%     ix = 13150; % location of interest
%     L(ix) = true;
%     D = dependencemap(M,L);
%     D = +D;
%     D(L) = D(L)+2;
%     imagesc(X(1,:),Y(:,2),D); axis image; axis xy
%     hold on
%     contour(X,Y,dem,30,'k')
%     caxis([0 2])
%
%
% See also: FLOWDIR, INFLUENCEMAP
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 15. March, 2009



% check input
siz = size(I);
nrc = numel(I);

if size(M,1) ~= nrc || size(M,2) ~= nrc
    error('TopoToolbox:incorrectinput',...
          'M or L have wrong dimensions')
end

I = I*1000;

I = (speye(nrc)-M)\I(:);
I = I>0;
I = reshape(I,siz);



