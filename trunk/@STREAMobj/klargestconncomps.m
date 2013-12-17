function S = klargestconncomps(S,k)

% retain k largest connected components in an instance of STREAMobj
%
% Syntax
%
%     S2 = klargestconncomps(S,k)
%
% Description
%
%     klargestconncomps extracts the k largest strongly connected 
%     components from the stream network in an instance of STREAMobj.
%     Largest refers to the number of edges (links between individual
%     pixels) in the network.
%
% Input arguments
%
%     S        instance of STREAMobj
%     k        number of largest components (default k=1)
%
% Output arguments
%
%     S2       new instance of STREAMobj that contains only the k largest
%              components of S
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 30. January, 2013


% check input arguments
narginchk(1,2)
if nargin == 1;
    k = 1;
else
    validateattributes(k,{'numeric'},{'scalar','integer','>',0},'klargestconncomps','k',2);
end

% create sparse adjacency matrix
nrc = numel(S.x);
M = sparse(double(S.ix),double(S.ixc),true,nrc,nrc);
M = M | M' | speye(nrc);
% use dmperm to find connected components
% see http://www.mathworks.de/matlabcentral/fileexchange/21366-findcomponents
% by Tim Davis
[~,p,~,r] = dmperm(M);
nc = length(r) - 1;

% label matrix
L = zeros(nrc,1);
% sort to find the largest conncompoments
[~,dd] = sort(diff(r),'descend');
if k>nc
    warning('TopoToolbox:STREAMobj',...
            ['There are only ' num2str(nc) ' connected components in the stream network']);
end
for tt = dd(1:min(nc,k)); %1:min(nc,k);
    L(p(r(tt):r(tt+1)-1)) = tt;
end

% adapt new STREAMobj to the reduced network
L = L>0;
I = L(S.ix);
S.ix = S.ix(I);
S.ixc = S.ixc(I);

IX = cumsum(L);
% IX([L(1); diff(IX)==0]) = 0;
S.ix  = IX(S.ix);
S.ixc = IX(S.ixc);

% II = S.ixc == 0;
% S.ix(II) = [];
% S.ixc(II) = [];

S.x = S.x(L);
S.y = S.y(L);
S.IXgrid = S.IXgrid(L);

