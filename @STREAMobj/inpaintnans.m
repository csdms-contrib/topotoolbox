function z = inpaintnans(S,DEM,extrap)

%INPAINTNANS inpaint missing values (nans) in a node attribute list
%
% Syntax
%
%     a = inpaintnans(S,A)
%     a = inpaintnans(S,an)
%     a = inpaintnans(...,extrap)
%
% Description
%
%     This function interpolates missing values of a variable measured
%     along a stream network using laplace's equation (see reference below).
%
% Input arguments
%
%     S      stream network (STREAMobj)
%     A      GRIDobj (must have the same coordinate system as S, but needs  
%            not be spatially aligned)
%     an     node attribute list with missing values
%     extrap true or false. 
%
% Output arguments
%  
%     a      node attribute list with missing values filled
%
% Example
%     
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','carve');
%     S   = STREAMobj(FD,'minarea',1000);
%     % now let's make a copy of DEM and corrupt it with many nans
%     DEM2 = DEM;
%     DEM2.Z(randperm(prod(DEM2.size),300000)) = nan;
%     % inpainting
%     z = inpaintnans(S,DEM2);
%     % and plotting
%     ax(1) = subplot(1,2,1);
%     plotdz(S,DEM2,'color','k')
%     ax(2) = subplot(1,2,2);
%     plotdz(S,z,'color','r')
%     linkaxes(ax,'xy');
%     % now zoom into the axes and move around to explore the results
%
%     
% See also: GRIDobj/inpaintnans
%
% Reference: 
% http://blogs.mathworks.com/steve/2015/06/17/region-filling-and-laplaces-equation/
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 17. August, 2015

% get node attribute list with elevation values
if isa(DEM,'GRIDobj')
    tf = validatealignment(S,DEM);
    if tf
        z = getnal(S,DEM);
    else
        z = interp(DEM,S.x,S.y);
    end
elseif isnal(S,DEM)
    z = DEM;
else
    error('Imcompatible format of second input argument')
end

if nargin == 2;
    extrap = false;
end

z = double(z(:));
inan = isnan(z);

if ~any(inan)
    return
end

% upstream distance
d  = S.distance;
% nr of nodes
nr = numel(S.IXgrid);

% create sparse adjacency matrices weighted by the inverse distance between
% nodes
A = sparse(double(S.ix(inan(S.ix))),double(S.ixc(inan(S.ix))),...
    1./(d(S.ix(inan(S.ix)))-d(S.ixc(inan(S.ix)))),nr,nr); 

A2 = sparse(double(S.ixc(inan(S.ixc))),double(S.ix(inan(S.ixc))),...
    1./(d(S.ix(inan(S.ixc)))-d(S.ixc(inan(S.ixc)))),nr,nr); 

% since there may be more than 1 upstream neighbor, A2 is normalized such
% that the sum of upstream weights equals 1-downstream weight. 
s = 1./(sum(spones(A2),2));
s(~inan) = 0;
A2 = spdiags(s,0,nr,nr)*A2;

% add A and A2
A = A+A2;
% and normalize again such that rows sum to one
s = 1./(sum(A,2));
s(~inan) = 0;
A = spdiags(s,0,nr,nr)*A;

% fill main diagonal with ones
A = speye(nr)-A;

% solve  
z(inan) = 0;
if ~extrap
    I = (streampoi(S,'channelhead','logical') | streampoi(S,'outlet','logical')) & inan;
    z(I) = nan;
end
z = A\z;

