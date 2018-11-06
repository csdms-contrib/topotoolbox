function [V,IX] = mapfromnal(FD,S,nal)

%MAPFROMNAL map values from node-attribute list to nearest upstream grid
%
% Syntax
%
%     V = mapfromnal(FD,S,nal)
%     [V,IX] = mapfromnal(FD,S,nal)
%
% Description
%
%     mapfromnal takes a STREAMobj S and an associated node-attribute list
%     nal and maps the values in the nal to the nearest grid values
%     measured along flowpaths based on the FLOWobj FD. S should be have
%     been derived from FD.
%
% Input arguments
%
%     FD      FLOWobj
%     S       STREAMobj
%     nal     node-attribute list
%
% Output arguments
%
%     V       GRIDobj with values derived from nal
%     IX      matrix with size V.size with linear indices into 
%             node-attributes of S. Elements in IX with no downstream
%             stream pixel are zero.
%
% See also: FLOWobj, STREAMobj, flowdistance, vertdistance2stream
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 6. November, 2018

narginchk(3,3)

% If the variable nal is a GRIDobj than extract the nal
if isa(nal,'GRIDobj')
    nal = getnal(S,nal);
end


nrnal = numel(S.x);
nalix = (1:nrnal)';

I   = STREAMobj2GRIDobj(S);
IX  = zeros(I.size);
IX(S.IXgrid) = nalix;

ix  = FD.ix;
ixc = FD.ixc;

for r = numel(ixc):-1:1
    if ~(I.Z(ixc(r)) && I.Z(ix(r)))
        IX(ix(r)) = IX(ixc(r));
    end
end

V   = GRIDobj(I,'single');
V.Z(:,:) = nan('single');
I   = IX~=0;
V.Z(I) = nal(IX(I));



