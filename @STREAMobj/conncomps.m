function [L,nc] = conncomps(S)

%CONNCOMPS labels of connected components (individual trees) in a stream network
%
% Syntax
%
%     [L,nc] = conncomps(S)
%
% Description
%
%     conncomps returns a node attribute list (nal) that contains labels
%     for connected components in a stream network. 
%
% Input arguments
%
%     S      STREAMobj
%
% Output arguments
%
%     L      node attribute list (nal) with labels
%     nc     number of connected components
%
% Example
%   
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     L = conncomps(S);
%     plotc(S,L)
%
% See also: STREAMobj, STREAMobj/klargestconncomps,
%           STREAMobj/extractconncomps
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 27. February, 2015 

nrc = numel(S.x);
M = sparse(double(S.ix),double(S.ixc),true,nrc,nrc);

[~,p,~,r] = dmperm(M | M' | speye(nrc));
nc = length(r) - 1;

% label matrix
L = zeros(nrc,1);
for tt = 1:nc;
    L(p(r(tt):r(tt+1)-1)) = tt;
end