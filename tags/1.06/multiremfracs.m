function M = multiremfracs(M,cutoff)

% remove drainage fractions in a multiple flow direction matrix
%
% Syntax
%
%     M = multiremfracs(M,cutoff)
% 
% Description
%
%     Sometimes the multiple flow direction matrix produces too dissipative
%     flow. This can be avoided by removing small fractions in the flow
%     direction matrix. Fractions are ratios of one transferred from each
%     cell to the neighboring cells. 
%
% Input
%
%     M         multiple flow direction matrix
%     cutoff    cutoff value (scalar). cutoff must have range ]0 .125[
% 
% Output
%   
%     Mr        reduced multiple flow direction matrix
%
% Example
%
%     load exampleDEM
%
%     M = flowdir(X,Y,dem);
%     L = false(size(dem));
%     ix = 10365; % location of interest
%     L(ix) = true;
%     I = influencemap(M,L);
%     I = +I;
%     I(L) = I(L)+2;
% 
%     subplot(1,2,1);
%     imagesc(X(1,:),Y(:,2),I); axis image; axis xy
%     hold on
%     contour(X,Y,dem,30,'k')
%     caxis([0 2])
%     title('influence map with multiple flow direction matrix')
% 
%     M = multiremfracs(M,0.12);
%     I = influencemap(M,L);
%     I = +I;
%     I(L) = I(L)+2;
%     subplot(1,2,2);
%     imagesc(X(1,:),Y(:,2),I); axis image; axis xy
%     hold on
%     contour(X,Y,dem,30,'k')
%     caxis([0 2])
%     title('influence map with reduced multiple flow direction matrix')
%
%
% References:
%
% Gruber, S., Huggel, C., Pike, R. 2009: Modelling mass movements and
% landslide susceptibility. In: Geomorphometry. Concepts, Software and
% Applications, Hengl, T., Reuter, H.I. (eds.), Elsevier, pp.527-550. 
%
%
% See also: FLOWDIR
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 15. March, 2009


% parse input data
p = inputParser; 
addRequired(p, 'cutoff', ...
   @(x)validateattributes(x, {'numeric'}, {'scalar','>',0,'<',0.125}, 'multiremfracs','cutoff',2)); 
parse(p,cutoff)


nrc = size(M,1);

II = max(M,[],2)<cutoff;

I = M>=cutoff;
M = (spdiags(II,0,nrc,nrc)+I).*M;

% Row standardization of M only necessary when multiple flow dir
M = spdiags(spfun(@(x) 1./x,sum(M,2)),0,nrc,nrc) * M;


