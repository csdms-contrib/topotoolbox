function S = trunk(S)

% extract trunk stream (longest stream) 
%
% Syntax
%
%     S2 = trunk(S)
%
% Description
%
%     trunk reduces a stream network to the longest streams in each stream
%     network tree (e.g. connected component). The algorithm identifies
%     the main trunk by sequently tracing the maximum downstream 
%     distance in upstream direction. 
%
% Input 
%
%     S    stream network (STREAMobj)
% 
% Output
%
%     S2   stream network (STREAMobj) with only trunk streams in each
%          connecoted component
%
% Example
%
%
% See also: chiplot, FLOWobj/flowpathextract, STREAMobj/klargestconncomps
%  
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 3. May, 2013

narginchk(1,1);
% downstream distance
nrc = numel(S.x);
dds = zeros(nrc,1);
% intercell distance
dx  = sqrt((S.x(S.ix)-S.x(S.ixc)).^2 + (S.y(S.ix)-S.y(S.ixc)).^2);
for r = 1:numel(S.ix);
    dds(S.ixc(r)) = max(dds(S.ixc(r)),dds(S.ix(r))+dx(r));
end

D = sparse(double(S.ix),double(S.ixc),dds(S.ix)+1,nrc,nrc);
OUTLET = any(D,1)' & ~any(D,2);
[~,Imax] = max(D,[],1);
II = false(nrc,1);
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


    