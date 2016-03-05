function [xy,Adj] = label2poly(L,X,Y,c)

% plot region outlines with polyline
%
% Syntax
%
%     label2poly(L)
%     label2poly(L,X,Y)
%     [xy,Adj] = label2poly(...)
%
% Description
%
%     label2poly extracts the outlines of regions in a label matrix and 
%     plots them. The outlines are drawn along the pixel edges. Without
%     output arguments, the function plots the outlines.
%
% Input arguments:
%
%     L         label matrix as returned by bwlabel or labelmatrix. 
%               A logical matrix can be supplied, too.
%     X,Y       coordinate matrices as returned by meshgrid
%
% Output arguments:
%
%     xy        coordinates
%     Adj       adjacency matrix as used for gplot
%
% Example:
%
%     [X,Y,P] = peaks(100);
%     L = bwlabel(P<-0.5);
%     subplot(1,2,1);
%     imagesc(X(1,:),Y(:,1),L);
%     axis image; axis xy
%     [xy,Adj] = label2poly(L,X,Y);
%     subplot(1,2,2);
%     gplot(Adj,xy,'k'); axis image
%
% see also: ixneighbors, bwlabel, labelmatrix, gplot
%
% Date: 3. August, 2011
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)



if nargin == 1 || isempty(X) || isempty(Y);
    [X,Y] = meshgrid(1:size(L,2),1:size(L,1));
end
if nargin<4;
    c = 'k';
end

% set nans in label matrix to zero
L(isnan(L)) = 0;

% general values needed later on
cs   = abs(Y(1)-Y(2));
cs2  = cs/2;
minX = min(X(1,:));
minY = min(Y(:,1));
maxX = max(X(1,:));
maxY = max(Y(:,1));

% identify border pixels
NHOOD = ones(3);
I = (imdilate(L,NHOOD) ~= L) | (imerode(L,NHOOD) ~= L);
L = double(L);

% set non-border pixel to nan
L(~I) = nan;

% call ixneighbor, which will connect all border pixels within an
% 8-neighborhood
[neighb(:,1) neighb(:,2)] = ixneighbors(L,I,8);

% remove double edges 
neighb = unique(sort(neighb,2),'rows');

% remove edges to neighbors with the same label
neighb(L(neighb(:,1))==L(neighb(:,2)),:) = [];

% calculate the coordinates of the labelgrid by averaging the coordinates
% along the neighbor edges
xy = [(X(neighb(:,1))+X(neighb(:,2)))/2  (Y(neighb(:,1))+Y(neighb(:,2)))/2];

% again remove double entries
xy = unique(xy,'rows');

% create coordinates of the label polygons along the boundaries of L 
xgr = (minX-cs2:cs2:maxX+cs2)'; 
nrx = numel(xgr);
ygr = (minY-cs2:cs2:maxY+cs2)'; 
nry = numel(ygr);
xy = [xy; ...
     [xgr repmat(minY-cs2,nrx,1)]; ...
     [xgr repmat(maxY+cs2,nrx,1)]; ...
     [repmat(minX-cs2,nry,1) ygr]; ...
     [repmat(maxX+cs2,nry,1) ygr]];

% number of label polygon coordinates
nrxy = size(xy,1);

% round coordinates to integer values to avoid floating point problems with ismember
ij  = [round((xy(:,1)-min(xy(:,1)))/cs2)+1 round((xy(:,2)-min(xy(:,2)))/cs2)+1];

% Construct node adjacency
[I1,loc1] = ismember(bsxfun(@minus,ij,[0 1]),ij,'rows');
[I2,loc2] = ismember(bsxfun(@minus,ij,[1 0]),ij,'rows');
i = [find(I1);find(I2)];
j = [loc1(I1);loc2(I2)];
Adj = sparse(i,j,1,nrxy,nrxy);

% if no output arguments, plot using gplot
if nargout == 0;
    label2polygplot(Adj,xy,gca,c);
    clear xy
elseif nargout == 1;
    xy = label2polygplot(Adj,xy);
end

end



function h = label2polygplot(Adj,xy,ax,c)


[i,j] = find(Adj);
[~, p] = sort(max(i,j));
i = i(p);
j = j(p);

X = [ xy(i,1) xy(j,1)]';
Y = [ xy(i,2) xy(j,2)]';
X = [X; NaN(size(i))'];
Y = [Y; NaN(size(i))'];

h = plot(ax,X(:),Y(:),[c '-']);

end
