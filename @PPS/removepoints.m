function P = removepoints(P,data)

%REMOVEPOINTS Remove points in point pattern
%
% Syntax
%
%     Pr = removepoints(P,tf)
%
% Description
%
%     removepoints removes points from a point pattern based on marks,
%     covariates or GRIDobjs. Points true in tf are removed.
%
% Input arguments
%
%     P       instance of PPS
%     data    logical node-attribute list, GRIDobj or vector same size as
%             npoints(P)x1
%
% Output arguments
%
%     Pr      instance of PPS with points removed.
%
%
% See also: PPS
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. February, 2019

if npoints(P) == numel(data)
    tf = data;
    P.PP(tf) = [];
else
    c = getcovariate(P,data);
    tf = any(c,2);
    P.PP(tf(P.PP)) = [];
end

