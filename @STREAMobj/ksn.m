function k = ksn(S,DEM,A,theta)

%KSN normalized steepness index
%
% Syntax
%
%     k = ksn(S,DEM,A,theta)
%     k = ksn(S,z,a,theta)
%
% Description
%
%     KSN returns the normalized steepness index. 
%
% Input arguments
%
%     S      STREAMobj
%     DEM    digital elevation model (GRIDobj)
%     A      flow accumulation as returned by flowacc (GRIDobj)
%     theta  concavity (e.g. 0.45)
%     z      node attribute list of elevation values
%     a      node attribute list of flow accumulation values
%
% Output arguments
%
%     k      normalized steepness index
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S);
%     DEM = imposemin(S,DEM);
%     A = flowacc(FD);
%     k = ksn(S,DEM,A,0.45);
%     k = smooth(S,k,'K',100);
%     subplot(2,1,1);
%     imageschs(DEM,DEM,'colormap',[.9 .9 .9],'colorbar',false);
%     hold on
%     plotc(S,k)
%     
%     subplot(2,1,2);
%     plotdz(S,DEM,'color',k);
%
% See also: STREAMobj/crs, STREAMobj/smooth, FLOWobj/flowacc
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 18. August, 2017

narginchk(3,4);
if nargin == 3
    theta = 0.45;
end

% get node attribute list with elevation values
if isa(DEM,'GRIDobj')
    validatealignment(S,DEM);
    z = double(getnal(S,DEM));
elseif isnal(S,DEM)
    z = double(DEM);
else
    error('Imcompatible format of second input argument')
end

% get node attribute list with flow accumulation values
if isa(A,'GRIDobj')
    validatealignment(S,A);
    a = double(getnal(S,A));
elseif isnal(S,DEM)
    a = double(DEM);
else
    error('Imcompatible format of second input argument')
end

% minima imposition to avoid negative gradients
z = imposemin(S,z);
% calculate gradient
g = gradient(S,z);
% upslope area
a = a.*S.cellsize.^2;

k = g./(a.^(-theta));




