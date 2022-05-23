function [v,nIX] = extractvaluesaroundpoints(P,z,varargin)

%EXTRACTVALUESAROUNDPOINTS Extract values around points 
%
% Syntax
%
%     v = extractvaluesaroundpoints(P,z)
%     v = extractvaluesaroundpoints(p,z,'pn',pv,...)
%     [v,nIX] = ...
%
% Description
%
%     extractvaluesaroundpoints extracts values from a GRIDobj or
%     node-attribute list z along the stream network around the points in
%     P. P is an instance of PPS and stores the stream network and points
%     located on the stream network. 
%
%     The function can be, for example, used to extract river gradients (or
%     ksn values) upstream and downstream of knickpoints, the difference of
%     which can be used as a metric of knickpoint prominence.
%
% Input arguments
%
%     P     instance of PPS
%     z     GRIDobj or node-attribute list
%     
%     Parameter name/value pairs
%
%     'direction'   {'both'}, 'upstream' or 'downstream'.
%     'dfrompoint'  distance from points
%     'aggfun'      anonymous aggregation function. {@mean}
%     'distance'    By default, distances between are calculated as the
%                   euclidean distance between stream nodes. Yet, other
%                   distance metrics (e.g. chi) can be used as well. Note
%                   that the value in dfrompoint must be given in the same
%                   distance units as those provided by here.
%    
% Output arguments
%
%     v       extracted value
%     nIX     cell array of linear indices into node-attribute list for all 
%             neighbors of all points in P.
%    
% Example
%
% See also: PPS, PPS/getmarks, nearest
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 21. May, 2022

p = inputParser;
p.FunctionName = 'extractvaluesaroundpoints';

addParameter(p,'direction','both',@(x) ischar(validatestring(x,{'both','upstream','downstream'})))
addParameter(p,'dfrompoint',10*P.S.cellsize)
addParameter(p,'distance',P.S.distance)
addParameter(p,'aggfun',@mean)

parse(p,varargin{:})

if isa(z,'GRIDobj')
    z = double(getnal(P.S,z));
else
    if isnal(P.S,z)
        z = double(z);
    else
        error('Cannot handle second input. It must be a GRIDobj or node-attribute list.')
    end
end

% get distance
d = p.Results.distance;
% and calculate inter-node distances as weights
w = abs(d(P.S.ix) - d(P.S.ixc));

direction = validatestring(p.Results.direction,{'both','upstream','downstream'});
switch direction
    case 'downstream'
        G = digraph(P.S.ix, P.S.ixc,w);
    case 'upstream'
        G = digraph(P.S.ixc, P.S.ix, w);
    case 'both'
        G = graph([P.S.ix; P.S.ixc],[P.S.ixc; P.S.ix], [w;w]);
end

% Calculate nearest neighbor indices using cellfun
pIX = num2cell(P.PP);

if nargout == 2
    nIX = cellfun(@(s) nearest(G,s,p.Results.dfrompoint),pIX,'UniformOutput',false);
    v   = cellfun(@(ix) p.Results.aggfun(z(ix)),nIX);
elseif nargout == 1
    v   = cellfun(@(s) p.Results.aggfun(z(nearest(G,s,p.Results.dfrompoint))),...
                pIX,'UniformOutput',true);
end

    
