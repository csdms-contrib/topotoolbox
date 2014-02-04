function G = STREAMobj2GRIDobj(S)

% convert an instance of STREAMobj to an instance of GRIDobj
%
% Syntax 
%
%     G = STREAMobj2GRIDobj(S)
%
% Description
%     
%     STREAMobj2GRIDobj converts an instance of STREAMobj S to a new instance
%     of GRIDobj G with the same spatial reference and extent as the instance
%     from which S has been derived.
%     
% Input arguments
%
%     S     STREAMobj
%     G     GRIDobj (contains logical raster with ones at stream locations
%           and zeros elsewhere.
%
%
% See also: STREAMobj, GRIDobj
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 22. March, 2013

narginchk(1,1)

G = GRIDobj([]);

props_G = fieldnames(G);
props_S = fieldnames(S);

[props_joint,b] = ismember(props_S,props_G); 

six = find(props_joint);
gix = b(six);

for r = 1:numel(six);
    G.(props_G{gix(r)}) = S.(props_S{six(r)});
end

G.Z = false(G.size);
G.Z(S.IXgrid) = true;

G.name = 'stream grid';





