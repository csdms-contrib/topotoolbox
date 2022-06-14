function [P,a,locb] = removeduplicates(P)

%REMOVEDUPLICATES removes duplicate points
%
% Syntax
%
%     P2 = removeduplicates(P)
%     [P2,a,locb] = ...
%
% Description
%
%     removeduplicates removes duplicate points in the point pattern P by
%     retaining only one of the duplicates.
%
% Input arguments
%
%     P      point pattern on stream network (class PPS)
%
% Output arguments
%
%     P      point pattern on stream network (class PPS)
%
% See also: PPS, PPS/hasduplicates
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. February, 2019

[u,a,b] = unique(P.PP,'stable');

if numel(u) == numel(P.PP)
    return
end

if nargout >= 2
    ix = (1:npoints(P))';
    ixp = accumarray(b,ix,[numel(u) 1],@(x) {x});
    ixp = cellfun(@(x) sort(x),ixp,'UniformOutput',false);
    locb = cellfun(@(x) x(2:end),ixp,'UniformOutput',false);
end

P.PP = P.PP(a);


