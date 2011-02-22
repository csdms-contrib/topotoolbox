function [I,Mstreams] = influencemap(M,I)

% downslope area for specific locations in a digital elevation model
%
% Syntax
%
%     [I,Mstreams] = influencemap(M,L)
%
% Description
%
%     influencemap returns a logical matrix masking the downslope part
%     of the digital elevation model that is drained by specified cells.
%
% Input
%
%     M         multiple or single flow direction matrix
%     L         logical matrix same size as the digital elevation from
%               which the flow direction matrix was derived from. true
%               values in L indicate cells whose downslope area is to be
%               calculated.
%
% Output
%
%     I         logical dependence matrix
%     Mstreams  modified flow direction matrix with nonzero entries only
%               for stream cells (may be used for plotting streams using
%               gplot (e.g. gplot(Mstreams,[X(:) Y(:)]);).  
%
% Example
%
%     load exampleDEM
%     M = flowdir(X,Y,dem,'exponent',5);
%     M = multiremfracs(M,0.05);
%     L = false(size(dem));
%     ix = 10360; % location of interest
%     L(ix) = true;
%     I = influencemap(M,L);
%     I = +I;
%     I(L) = I(L)+2;
%     imagesc(X(1,:),Y(:,2),I); axis image; axis xy
%     hold on
%     contour(X,Y,dem,30,'k')
%     caxis([0 2])
%
%
% See also: FLOWDIR, DEPENDENCEMAP
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 16. September, 2009





% check input
siz = size(I);
nrc = numel(I);

if size(M,1) ~= nrc || size(M,2) ~= nrc
    error('TopoToolbox:incorrectinput',...
          'M or L have wrong dimensions')
end

I = I*1000;

I = (speye(nrc)-M')\I(:);
I = I>0;

if nargout == 2
    Mstreams = spdiags(I,0,nrc,nrc)*M;
end

I = reshape(I,siz);



