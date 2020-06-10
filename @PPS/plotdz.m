function varargout = plotdz(P,varargin)

%PLOTDZ plot upstream distance version elevation or covariate of a PPS
%
% Syntax
%
%     plotdz(P)
%
% Description
%
%     plot distance versus elevation of points and stream network in an
%     instance of PPS.
%
% Input arguments
%
%     P       instance of PPS (needs z-property unless parameter z is set
%             to different node-attribute list)
%
%     Parameter name/value pairs
%
%     all parameters that work for STREAMobj/plotdz
%
%     all parameters that work for scatter
%
%     'z'     default is the node-attribute of elevation values stored in
%             the PPS object. Takes any other node-attribute list.
%
%
% See also: PPS, PPS/npoints 
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. February, 2019

p = inputParser;
p.KeepUnmatched = true;
addParameter(p,'z','z');
addParameter(p,'Marker','o');
addParameter(p,'MarkerEdgeColor','k');
addParameter(p,'MarkerFaceColor','w');
addParameter(p,'SizeData',20);
addParameter(p,'MarkerFaceAlpha',0.5);
addParameter(p,'distance',[],@(x) isnal(P.S,x));

% Parse
parse(p,varargin{:});

UnMatched = p.Unmatched;
UnMatched.distance = p.Results.distance;

tf = ishold;

zz  = getcovariate(P,p.Results.z);

hl = plotdz(P.S,zz,UnMatched);
hold on
if isempty(p.Results.distance)
    d  = P.S.distance;
else
    d  = p.Results.distance;
end
Results = p.Results;
Results = rmfield(Results,{'distance' 'z'});

% Does the field 'dunit' exist
if isfield(UnMatched,'dunit')
    switch UnMatched.dunit
        case 'km'
            d = d/1000;
    end
end
pnpv    = expandstruct(Results);

hs = scatter(d(P.PP),zz(P.PP),pnpv{:});
if ~tf
    hold off
end

if nargout >= 1
    varargout{1} = hl;
end
if nargout == 2
    varargout{2} = hs;
end
end

function pnpv = expandstruct(s)

pn = fieldnames(s);
pv = struct2cell(s);

pnpv = [pn pv];
pnpv = pnpv';
pnpv = pnpv(:)';
end
