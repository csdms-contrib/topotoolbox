function H = sbplot(S,X,Y,v)

% plot subbasins and connectivity between subbasin gauges
%
% Syntax
%
%     sbplot(S,X,Y)
%     sbplot(S,X,Y,v)
%     h = sbplot(...)
%
% Description
%
%     sbplot visualizes sub-basins and their connectivity
%
% Input
%
%     S     sub-basin structure array as returned by sbstruct
%     X,Y   coordinate matrices of the DEM
%     v     vector with values for each subbasin
%
% Output
%
% Example
%
%     load exampleDEM
%     M = flowdir(X,Y,dem,'type','single');
%     A = flowacc(M,size(dem));
%     W = A>100;
%     [nodes,nodes] = streamorder(M,W);
%     S = sbstruct(M,nodes);
%     P = sbprops(S,X,Y,false);
%     sbplot(S,X,Y,[P.Area])
%
%
% See also: COORD2IND, SBPLOT, SBPROPS
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 30. October, 2009


error(nargchk(3, 4, nargin))
if nargin == 3;
    plotvar = false;
else
    plotvar = true;
    if numel(v)~=numel(S.SBIx)
        error('TopoToolbox:incorrectinput',...
          'v must have same size as S.SBIx')
    end
end

siz   = size(X);

% nr of sub-basins
nSBIx = numel(S.SBIx);

% create label matrix
L   = nan(siz);
for r = 1:nSBIx;
    idx = S.SBIdxList{r};
    if plotvar
        L(idx) = v(r);
    else
        L(idx) = r;
    end
end

h = imagesc(X(1,:),Y(:,1),L);
axis image
axis xy

hold on
gplot(S.Adj,[X(S.OutletIdx) Y(S.OutletIdx)],'k-');
plot(X(S.OutletIdx),Y(S.OutletIdx),'k.')

% text(X(S.OutletIdx),Y(S.OutletIdx),...
%      cellfun(@(x) num2str(x),num2cell(S.SBIx),'UniformOutput',false),...
%      'FontSize',14,...
%      'FontWeight','b');

 
hold off

if nargout == 0;
else
   H = h;
end


