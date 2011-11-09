function [U,R,c] = upslopesidelength(M,I,cs)

% calculate the upslope side length and connectivity
%
% Syntax
%
%     [U,u,c] = upslopesidelength(M,I,cs)
%
% Description
%
%     The upslope side length of a vegetated patch is a useful indicator of
%     the patch's capacity to collect water that flows downslope (Imeson &
%     Prinsen, 2004).
%
% Input
%
%     M         sparse flow direction matrix (single and multi)
%     I         vegetation patches/vegetation density
%     cs        cell size
%
% Output
% 
%     U         Pixel with upslope bare areas (vegetated pixels that 
%               receive more than a fraction of 0.5 from surrounding bare
%               areas)
%     u         mean vegetated patch upslope side length
%     c         connectivity between bare areas (fraction)
%
% Example
%
%     load exampleDEM
%     M = flowdir(X,Y,dem,'exponent',5);
%     A = flowacc(M,size(dem));
%     surf(X,Y,dem,log(A)) 
%
% References
%     
%     Imeson, A.C., Prinsen, H.A.M, 2004. Vegetation patterns as biological
%     indicators for identifying runoff and sediment source and sink areas
%     for semi-arid landscapes in Spain. Agriculture, Ecosystems and
%     Environment, 104, 333-342.
% 
%
% See also: FLOWDIR, CROSSFLATS, UPSLOPESTATS
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 15. March, 2009


siz = size(I);
I(isnan(I)) = 0;

% calculate the total number of upslope side cells
U = (M'* (~I(:))) > 0.4 & (I(:)>0); 
U = reshape(U,siz);
% total number of patches
CC = bwconncomp(I>0);
numpatches = CC.NumObjects;

% upslope side length of each patch
L = zeros(siz);
R = zeros(siz);
for r = 1:numpatches;
    L(CC.PixelIdxList{r}) = sum(U(CC.PixelIdxList{r}))*cs;
    R(CC.PixelIdxList{r}) = sum(U(CC.PixelIdxList{r}))*cs/numel(CC.PixelIdxList{r})*cs.^2;
end

% mean upslope side length
u = sum(U(:))*cs/numpatches;

% connectivity
c = flowacc(M,[],1-I)./flowacc(M,size(I));
c = mean(c(:));
