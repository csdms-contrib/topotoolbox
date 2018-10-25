function z = cumtrapz(S,G)

%CUMTRAPZ Cumulative trapezoidal numerical integration along a stream network
%
% Syntax
%
%     z = cumtrapz(S,G)
%     z = cumtrapz(S,g)
%
% Description
%
%     cumtrapz(S,G) computes the cumulative integral of G with respect to
%     the distance along the stream network S using trapezoidal
%     integration. S must be a stream network. The second input must either
%     be a GRIDobj G or a node attribute list g.
%     
% Input arguments
%
%     S     STREAMobj
%     G     GRIDobj
%     g     node attribute list
%
% Output arguments
%
%     z     node attribute list
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     A = flowacc(FD);
%     a = getnal(S,A)*DEM.cellsize^2;
%     ghat = 1./(a.^0.45);
%     z = cumtrapz(S,ghat);
%     plotdz(S,z)
%
% See also: STREAMobj, STREAMobj/gradient
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 29. December, 2015


% get node attribute list with elevation values
if isa(G,'GRIDobj')
    validatealignment(S,G);
    g = getnal(S,G);
elseif isnal(S,G)
    g = G;
else
    error('Imcompatible format of second input argument')
end

z = zeros(size(g));
d = distance(S,'node_to_node');

ix  = S.ix;
ixc = S.ixc;

for r = numel(ix):-1:1
    z(ix(r)) = z(ixc(r)) + (g(ixc(r))+(g(ix(r))-g(ixc(r)))/2)*d(ix(r));
end

%% Is cumtrapz really working correctly. Here is the test:
% d  = S.distance;
% d2 = cumtrapz(S,ones(size(S.x)));
% tf = isequal(d,d2)


