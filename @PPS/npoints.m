function np = npoints(P)

%NPOINTS number of points in the point pattern
%
% Syntax
%
%     np = npoints(P)
%
% Description
%
%     npoints returns the number of points in the the point pattern.
%
% Input arguments
%
%     P      point pattern on stream network (class PPS)
%
% Output arguments
%
%     np     number of points (integer scalar)
%
%
% See also: PPS, PPS/tlength 
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. February, 2019


np = numel(P.PP);