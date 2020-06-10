function d = netdist(P,varargin)

%NETDIST Shortest network distance
%
% Syntax
%
%     d = netdist(P)
%
% Description
%
%     netdist computes the distance along a stream network. It calculates
%     the distance of each node in the network to the nearest point in the
%     PPS object P.
%
% Input arguments
%
%     P       PPS
% 
%     Parameter name/value pairs
%
%     'dir'     'both' (default) calculates the distance in upstream and
%               downstream direction. 'up' calculates distances only in
%               upstream direction and 'down' in downstream direction.
% 
%     'split'   false (default) or true. If true, netdist will do the
%               computation in parallel on each drainage basin. If the
%               Parallel Computing Toolbox is available, that may be 
%               faster. However, splitting the network in its connected
%               components produces some considerable computational
%               overhead.
%
% Output arguments
%
%     d     node-attribute list with distances. Nodes that cannot be
%           reached will have a distance of inf. 
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     S   = STREAMobj(FD,'minarea',500);
%     S   = klargestconncomps(S);
%     P   = PPS(S,'runif',100,'z',DEM);
%     d   = netdist(P);
%     plotc(P.S,d);
%
% See also: STREAMobj/netdist, PPS/voronoi
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 23. May, 2019


if npoints(P) == 0
    d = inf(size(P.S.x));
    return
end

d = netdist(P.S,P.S.IXgrid(P.PP),varargin{:});