function S = trunk(S,varargin)

%TRUNK extract trunk stream (longest stream) 
%
% Syntax
%
%     S2 = trunk(S)
%     S2 = trunk(S,A)
%     S2 = trunk(S,a)
%
% Description
%
%     trunk reduces a stream network to the longest streams in each stream
%     network tree (e.g. connected component). The algorithm identifies
%     the main trunk by sequently tracing the maximum downstream 
%     distance in upstream direction. 
%
%     If the second input argument is a flow accumulation grid A (or node
%     attribute list a) than the function will extract the main trunk by
%     tracing the maximum upstream area in upstream direction.
%
% Input 
%
%     S    stream network (STREAMobj)
% 
% Output
%
%     S2   stream network (STREAMobj) with only trunk streams in each
%          connecoted component
%     A    GRIDobj with flow accumulation values (as returned by the
%          function flowacc). 
%     a    node-attribute list (nal) (e.g. getnal(S,flowacc(FD))).
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     St = trunk(S);
%     plot(S)
%     hold on
%     plot(St)
%     legend('Stream network','Trunk rivers')
%
%
% See also: chiplot, FLOWobj/flowpathextract, STREAMobj/klargestconncomps
%  
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 8. September, 2019


narginchk(1,2);

if nargin == 2
    if isa(varargin{1},'GRIDobj')
        a = getnal(S,varargin{1});
    elseif isnal(S,varargin{1})
        a = varargin{1};
    else
        error('Cannot handle second input argument.');
    end
end
        
        

% downstream distance
nrc = numel(S.x);

if nargin == 1
    dds = distance(S,'max_from_ch');
else
    dds = a;
end

D        = sparse(double(S.ix),double(S.ixc),dds(S.ix)+1,nrc,nrc);
OUTLET   = any(D,1)' & ~any(D,2);
[~,Imax] = max(D,[],1);
II       = false(nrc,1);
II(Imax) = true;

I = false(nrc,1);
I(OUTLET) = true;
for r = numel(S.ix):-1:1
    I(S.ix(r)) = I(S.ixc(r)) && II(S.ix(r));
end

L = I;
I = L(S.ixc) & L(S.ix);

S.ix  = S.ix(I);
S.ixc = S.ixc(I);

IX    = cumsum(L);
S.ix  = IX(S.ix);
S.ixc = IX(S.ixc);

S.x   = S.x(L);
S.y   = S.y(L);
S.IXgrid   = S.IXgrid(L);
