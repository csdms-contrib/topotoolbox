function [S,ixnodes,M] = streamorder(M,W)

% calculate Strahler Stream Order from flow direction matrix
%
% Syntax
%
%     [S,nodes] = streamorder(M,W)
%     [S,nodes,Mstreams] = streamorder(M,W)
%
% Description
%
%     The Strahler Stream Order is a way to classify rivers based on a
%     hierarchy of tributaries. First-order streams don't have tributaries.
%     When two first-order streams come together, they form a second-order
%     stream. When two second-order streams come together, they form a
%     third-order stream. When two streams of different order confluence, 
%     they form a stream of maximum order of both.
%
%     streamorder returns the Strahler Stream Order based on a single flow 
%     direction matrix (M) and channel matrix (W). W is a logical matrix
%     and must have the same size as the digital elevation model from which 
%     the flow direction matrix has been calculated. It may either contain 
%     only channel starts or the the channel network.
% 
%     The output matrix S has the same size as W and contains the Strahler
%     Order for each cell. Non-channel cells are set to zero. nodes 
%     contains the linear indices of channel confluences. The third output 
%     is a sparse adjacency matrix which can be used to plot the stream 
%     network using gplot (gplot(Ms,[X(:) Y(:)], where X Y are the 
%     coordinate matrices for the Digital Elevation Model created by 
%     meshgrid).
%
% Input
%   
%     M         single flow direction raster
%     W         channel matrix (logical matrix, true where channels/
%               channelheads are) 
%
% Output
%
%     S         stream order matrix
%     nodes     linear index of cells where streams confluence
%     Mstreams  modified flow direction matrix with nonzero entries only
%               for stream cells (may be used for plotting streams using
%               gplot (e.g. gplot(Mstreams,[X(:) Y(:)]);).
%
% Example:
%
%     load exampleDEM
%     % calculate flow accumulation and direction
%     [A,M] = ezflowacc(X,Y,dem,'type','single');
%     % let's simply assume that channels start where
%     % A is larger than 100;
%     W = A>100;
%     % and calculate the strahler stream order ...
%     [S,nodes] = streamorder(M,W);
%     % ... and visualize it
%     subplot(1,2,1); 
%     pcolor(X,Y,+W); axis image; shading flat;
%     colorbar
%     title('Stream Network')
%     subplot(1,2,2);
%     pcolor(X,Y,S); axis image; shading flat;
%     colorbar
%     hold on
%     plot(X(nodes),Y(nodes),'ks','MarkerFaceColor','g')
%     title('Strahler Stream Order')
% 
% See also: EZFLOWACC, FLOWDIR_SINGLE
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 15. March, 2009




% get grid size
nrc = numel(W);
siz = size(W);

% The number of elements in W must be the same
% as the number of rows/cols in M
if nrc~= size(M,1) || nrc ~= size(M,2)
    error('There is a mismatch between M and W')
end

% if no channelheads/network are supplied
if ~any(W(:))
    S = zeros(siz);
    ixnodes = [];
    return
end
    
% check if single flow direction Matrix
% is used
if any(sum(spones(M),2)>1);
    error('single flow direction matrix must be used')
end

% force column vector
W = +W(:);

% detect, were channels are
% (only applies to when only channelstarts are given
% in matrix W. Since it is hard to guess if channelstarts
% or channel networks are supplied, the next two lines 
% are executed so or so).
T = (speye(nrc)-M')\W;
T = +(T~=0);

% remove values in M were no channels are
M = spdiags(T,0,nrc,nrc)*M;

% find channel heads
heads = ((sum(M,1) == 0)' & T);
ixheads = find(heads);
heads = +heads;
heads(ixheads) = ixheads;

% find channel nodes
nodes = (sum(M,1))';
nodes = max(nodes-1,0);
% nr of nodes
nrn = sum(nodes);
ixnodes = find(nodes);
nodes(ixnodes) = ixnodes;

% create sparse matrix to connect channel heads with nodes
[ic,icd,val]  = find(M);
I  = logical(nodes(icd));
val(I) = 0;
icnodes  = ic(I);
icdnodes = (1:numel(icnodes))'+nrc;
valnodes = ones(size(icdnodes));

% nr of rows and cols of distribution matrix
n  = nrc+nrn*2;

% D is a slightly different matrix than M. Some extra rows and 
% columns prevent information to be carried throughout the whole
% matrix. Instead, information on the indices is carried only to 
% respective nodes.
D  = sparse([ic;icnodes],[icd;icdnodes],[val;valnodes],n,n);

% Solve the equation and supply indices of channel heads and nodes
% on the right hand side
CONN = (speye(n)-D')\[heads+nodes;zeros(nrn*2,1)];

% IX2 (IX1) contain the linear indices of the upper (lower) nodes 
[IX2,IX1,IX1]  = find(CONN);

% The extra rows added to D now contain the information between
% nodes 
II = (IX2>nrc);
IX2(II) = icd(I);

% channel head or node cells IX1c drain to nodes IX2c
IX2c = IX2(II);
IX1c = IX1(II);

% if there are only channels of order 1
if isempty(IX2c)
    S = zeros(siz);
    S([IX1;IX2]) = 1;
    return
end

% channel heads and nodes that drain to their lower tributary partners
IX2(II) = [];
IX1(II) = [];

% now examine channel head/node to node connections
% create new linear indices of channels and 
% keep connectivity
[ixx,IX,IX] = unique([IX1c;IX2c]);

% new linear index for connections
nIX = (1:numel(ixx))';
IX  = reshape(IX,[],2);
nIX = nIX(IX);

% number of node connections
nrcc = max(nIX(:));

% each channel cell has value one at the beginning
W1  = ones(nrcc,1);

% set abort criterion 
undone = true;

% here follows the central function to determine stream order for
% links between nodes. It's a kind of recursive algorithm
while undone 
    % accumulate stream order using the function @strahler
    W2 = accumarray(nIX(:,2),W1(nIX(:,1)),[nrcc 1],@strahler,1);
    
    % the abortion criterion is set when there are no changes to
    % W2 anymore
    if W1==W2;
        undone = false;
    end   
    W1 = W2;
end

% now map the values back to the original extent
% of the dataset
S = zeros(siz);
S(IX2c) = W2(nIX(:,2));
S(reshape(W~=0,siz) & S==0) = 1;
S(IX2) = S(IX1);


% finito


function y = strahler(x)
% function for the recursive determination of the 
% strahler order for a node
if numel(x) == 1;
    y = x;
else
    if numel(unique(x))==1;
        y = x(1)+1;
    else
        y = max(x);
    end
end
    
    
    
    