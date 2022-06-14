function l = tlength(P)

%TLENGTH total length of the stream network
%
% Syntax
%
%     l = tlength(P)
%
% Description
%
%     tlength returns the total length of the stream network measured in
%     units of distance.
%
% 
% See also: PPS, PPS/npoints 
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. February, 2019

l = info(P.S,'totallength');