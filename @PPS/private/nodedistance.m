function d = nodedistance(P,varargin)

%NODEDISTANCE Compute edge list with node distances of the stream net
%
% Syntax
%
%     d = nodedistance(P)
%     d = nodedistance(P,pn,pv,...)
%
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 1. June, 2019

p = inputParser;
addParameter(p,'marks',[]);
addParameter(p,'d3d',false)
addParameter(p,'val',[]);

% Parse
parse(p,varargin{:});

if isempty(p.Results.val)
    if ~p.Results.d3d
        d = sqrt((P.S.x(P.S.ix)-P.S.x(P.S.ixc)).^2 + ...
            (P.S.y(P.S.ix)-P.S.y(P.S.ixc)).^2);
    else
        z = getcovariate(P,'z');
        d = sqrt((P.S.x(P.S.ix)-P.S.x(P.S.ixc)).^2 + ...
            (P.S.y(P.S.ix)-P.S.y(P.S.ixc)).^2 + ...
            (z(P.S.ix)-z(P.S.ixc)).^2);
    end
else
    c = getcovariate(P,p.Results.val);
    d = abs(c(P.S.ix) - c(P.S.ixc));
end
    
