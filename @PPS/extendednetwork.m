function [G,p,pback] = extendednetwork(P,varargin)
%EXTENTEDNETWORK Extend network to account for duplicate points
%
% Syntax
%
%     [G,p,pback] = extendednetwork(P)
%
% Description
%
%     extendednetwork creates a modified graph of the stream network stored
%     in a PPS object so that there are no duplicate points. 
%
% Input arguments
%
%     P      instance of PPS
%
% Output arguments
%
%     G      instance of graph
%     p      index of points in G
%     pback  index to map points back to P
%
%
% See also: PPS
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. February, 2019

dummydist = P.S.cellsize/10;

d = nodedistance(P);

if ~hasduplicates(P)
    G = graph(S.ix,S.ixc,d);
    return
end

[Pun,a,locb] = removeduplicates(P);

% count duplicates
ndup = sum(cellfun(@length,locb));

% common node

G = graph(P.S.ix,P.S.ixc,d);
n = numnodes(G);
G = addnode(G,ndup);

neigh = Pun.PP;
stw = cellfun(@(s,t) [s(:) repmat(t,numel(s),1) repmat(dummydist,numel(s),1)],...
              locb,num2cell(neigh),'UniformOutput',false);
stw = vertcat(stw{:});
s   = n+1:n+ndup;
s   = s(:);

G = addedge(G,s,stw(:,2),stw(:,3));

p = [Pun.PP;s];
pback = [a;stw(:,1)];





