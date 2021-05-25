function S = binarize(S,a)
%BINARIZE Make a stream network a strictly binary tree
%
% Syntax
%
%     S2 = binarize(S)
%     S2 = binarize(S,a)
%
% Description
%
%     Stream networks are commonly binary trees. This means that, when
%     moving upstream, a stream branches into a two streams. However, there
%     are cases where the stream network branches into three or more
%     streams which might be unwanted in some analysis.
%
%     The function binarize removes tributaries in confluences so that
%     there is a maximum of two tributaries. Tributaries to be removed are
%     the ones that have the smallest downstream distance, or those having
%     the smallest upstream area supplied as GRIDobj or node attribute list
%     a.
%
% Input arguments
%
%     S      STREAMobj
%     a      upstream area (GRIDobj or node-attribute list)
%
% Output arguments
%
%     S2     STREAMobj
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM);
%     S  = STREAMobj(FD,'minarea',50);
%     S2 = binarize(S);
%     plot(S,'b')
%     hold on
%     plot(S2,'r')
%     hold off
%
% See also: STREAMobj/streamorder, STREAMobj/modify
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 22. May, 2021
    

if nargin == 1
    a = distance(S,'max_from_ch');
else
    if isa(a,'GRIDobj')
        a = getnal(S,a);
    else
        if ~isnal(S,a)
            error('TopoToolbox:binarize',...
                'The second input argument must be a valid node-attribute list.')
        end
    end
end

% Calculate the number of upstream pixels
n = numel(S.x);

ix = accumarray(S.ixc,S.ix,[n 1],@(ix)getnonbinaryneighbors(ix));
ix = vertcat(ix{:});

I  = streampoi(S,'outl','logical');
II  = ismember(S.ix,ix);
S.ix(II) = [];
S.ixc(II) = [];


for r = numel(S.ix):-1:1
    I(S.ix(r)) = I(S.ixc(r));
end

S = subgraph(S,I);

% M = sparse(S.ix,S.ixc,1,n,n);
% nupstreamneighbors = sum(M,2);
% I = nupstreamneighbors>=3;

function ix = getnonbinaryneighbors(ix)

if numel(ix) <= 2
    ix = {[]};
    return
end

[~,ixx] = sort(a(ix),'descend');
ix = {ix(ixx(3:end))};
end

end



