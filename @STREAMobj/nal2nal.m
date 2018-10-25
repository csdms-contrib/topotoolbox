function nal2 = nal2nal(S2,S1,nal1,fillval)

%SETNAL 
%
% Syntax
%
%     nalB = nal2nal(SB,SA,nalA)
%     nalB = nal2nal(SB,SA,nalA,fillval)
%     nalB = nal2nal(SB,SA,nalA,nalB)
%
% Description
%
%     nal2nal maps the node-attribute list (nal) nalA of the stream network 
%     SA to another stream network SB. SA must be a subgraph of SB, i.e.,
%     SA is formed from a subset of the vertices (nodes) in SB.
%
% Input arguments
%
%     SB       STREAMobj
%     SA       STREAMobj that is a subgraph of SB
%     nalA     node-attribute list of SA
%     fillval  value to be assigned to nodes in nalB that are not members
%              of the network SA. By default, the value is NaN.
%     nalB     nal of SB. nal2nal will overwrite values at vertex-locations
%              in SA.
%     
% Output arguments
%
%     nalB     nal of SB
%
%
% See also: STREAMobj2GRIDobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 5. September, 2018


if nargin == 3
    fillval = nan;
end

nal2size = size(S2.IXgrid);

if isscalar(fillval)
    nal2 = repmat(fillval,nal2size);
else
    tf = isnal(S2,fillval);
    if ~tf
        error('TopoToolbox:setnal',...
            ['The forth input argument must be a scalar or\n'...
             'a node-attribute list of S2']);
    end
    nal2 = fillval;
end


[I,locb] = ismember(S1.IXgrid,S2.IXgrid);

if ~all(I)
    error('TopoToolbox:setnal',...
            ['S1 must be a subgraph of S2. All nodes in S1\n'...
             'must be a member of the network S2']);
end
    
nal2(locb) = nal1;