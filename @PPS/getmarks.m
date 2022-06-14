function m = getmarks(P,cova)
%GETMARKS Extract point marks
%
% Syntax
%
%     m = getmarks(P,data)
%
% Description
%
%     getmarks extracts values for each point in the point pattern P from
%     GRIDobjs or node-attribute lists cova.
%
% Input arguments
%
%     P       instance of PPS
%     data    cell array
%
% Output arguments
%
%     m       marks (vector or matrix with a row for each point in P)
%
%
% See also: PPS
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. February, 2019


c = getcovariate(P,cova);

if istable(c)
    
else
    m = c(P.PP,:);
end