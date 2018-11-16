function k = curvature(S,varargin)

%CURVATURE curvature or 2nd derivative of a STREAMobj
%
% Syntax
%
%     k = curvature(S)
%     k = curvature(S,DEM)
%     k = curvature(S,nal)
%     k = curvature(S,x,y)
%
% Description
%
%     curvature returns a node attribute list with curvature values. These
%     values can either represent the planform curvature (curvature in the
%     plane), profile curvature (second derivative of the longitudinal
%     stream profile) or second derivative of any other variable provided
%     as GRIDobj or node attribute list.
%
%     There is no straightforward way to calculate curvature at river
%     confluences since there are two or more pixels upstream. The 
%     algorithm uses the upstream pixel with the longest downstream flow
%     distance to derive curvature.
%
% Input arguments
%
%     S      STREAMobj
%     DEM    Digital elevation model (GRIDobj)
%     nal    node attribute list
%     x,y    node attribute list of coordinates 
%            curvature(S) is the same as curvature(S,S.x,S.y)
%
% Output arguments
%
%     k      node attribute list with curvature values
%
% Example
%
%     % calculate the curvature of a smoothed planform river network
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S  = STREAMobj(FD,'minarea',1000);
%     S  = klargestconncomps(S);
%     y = smooth(S,S.y,'k',1000,'nstribs',true);
%     x = smooth(S,S.x,'k',1000,'nstribs',true);
%     k = curvature(S,x,y);
%     plotc(S,k);
%     caxis(repmat(min(abs(caxis)),1,2).*[-1 1])
%     h = colorbar;
%     h.Label.String = 'Curvature';
%
%
% See also: STREAMobj/gradient, STREAMobj/smooth
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 25. January, 2016

narginchk(1,3)

% upstream distance
d = S.distance;

% This matrix contains the finite central differences
[I,loc] = ismember(S.ixc,S.ix);

%         i-1           i        i+1
colix  = [S.ixc(loc(I)) S.ixc(I) S.ix(I)];
rowix  = S.ixc(I);

val    = [2./((d(colix(:,2))-d(colix(:,1))).*(d(colix(:,3))-d(colix(:,1)))) ...
    -2./((d(colix(:,3))-d(colix(:,2))).*(d(colix(:,2))-d(colix(:,1)))) ...
    2./((d(colix(:,3))-d(colix(:,2))).*(d(colix(:,3))-d(colix(:,1))))];


dd = distance(S,'max_from_ch');
I = (dd(colix(:,2)) - dd(colix(:,3)))>=(sqrt(2*S.cellsize.^2)+S.cellsize/2);
%     colix(I,:) = [];
val(I,:) = 0;
%     rowix(I,:) = 0;


% second-derivative matrix
%     nrrows = size(colix,1);
rowix  = repmat(rowix,1,3);
nr     = numel(S.x);
Asd    = sparse(rowix(:),colix(:),val(:),nr,nr);

% check input arguments
if nargin == 1 || nargin == 3
    
    if nargin == 1
        x = S.x;
        y = S.y;
    else
        x = varargin{1};
        y = varargin{2};
    end
    
    dxx = Asd*x;
    dyy = Asd*y;
    
    dx  = gradient(S,x);
    dy  = gradient(S,y);
    
    % curvature in plane
    k   = (dx.*dyy - dy.*dxx)./(dx.^2 + dy.^2).^(3/2);
    k(isnan(k)) = 0;
    
else
    % get node attribute list with elevation values
    DEM = varargin{1};
    if isa(DEM,'GRIDobj')
        validatealignment(S,DEM);
        z = getnal(S,DEM);
    elseif isnal(S,DEM)
        z = DEM;
    else
        error('Imcompatible format of second input argument')
    end
    
    k      = Asd*z;
    k(isnan(k)) = 0;
end

