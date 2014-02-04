function [ix,ixc,frac] = find(FD)

% find indices and values of edges in the flow direction graph
%
% Syntax
%
%     [ix,ixc] = find(FD); % for single flow direction
%     [ix,ixc,frac] = find(FD); % for multiple flow direction
%
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 23. February, 2013


ix = FD.ix;
ixc = FD.ixc;
frac = FD.fraction;