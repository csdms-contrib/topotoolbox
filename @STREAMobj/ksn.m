function k = ksn(S,DEM,A,varargin)

%KSN normalized steepness index
%
% Syntax
%
%     k = ksn(S,DEM,A)
%     k = ksn(S,z,a)
%     k = ksn(...,theta)
%     k = ksn(...,theta,K)
%
% Description
%
%     KSN returns the normalized steepness index using a default concavity
%     index of 0.45.
%
% Input arguments
%
%     S      STREAMobj
%     DEM    digital elevation model (GRIDobj)
%     A      flow accumulation as returned by flowacc (GRIDobj). Note that
%            flowacc returns the number of pixels. The function ksn
%            calculates area in m^2 from these values internally.
%     z      node attribute list of elevation values
%     a      node attribute list of flow accumulation values
%     theta  concavity (default 0.45)
%     K      smoothing factor K (by default 0 = no smoothing). See function
%            STREAMobj/smooth for details
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
%     k = ksn(S,DEM,A,0.45,100);
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
% Date: 26. August, 2021


p = inputParser;
addOptional(p,'theta',0.45,...
    @(x) validateattributes(x,{'double','single'},{'scalar','positive'}));
addOptional(p,'smooth',0,...
    @(x) validateattributes(x,{'double','single'},{'scalar','positive'}));
parse(p,varargin{:});


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
elseif isnal(S,A)
    a = double(A);
else
    error('Imcompatible format of second input argument')
end

% minima imposition to avoid negative gradients
z = imposemin(S,z,0.00001);
% calculate gradient
g = gradient(S,z);
% upslope area
% if ~isgeographic(S)
%     a = a.*S.cellsize.^2;
% end

k = g./(a.^(-p.Results.theta));

if p.Results.smooth ~= 0
    k(k==0) = 0.0001;
    ks = smooth(S,log(k),'K',p.Results.smooth);
%     sig = var(ks-log(k));
    k  = exp(ks)*exp(var(-ks/2));
end


