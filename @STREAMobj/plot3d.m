function h = plot3d(S,DEM)

%PLOT3D 3D plot of a stream network
%
% Syntax
%
%     plot3d(S,DEM)
%     h = plot3d(S,DEM)
%
% Description
%
%     This function plots a 3D view of the stream network.
%
% Input arguments
%
%     S     stream network (STREAMobj)
%     DEM   digital elevation model (GRIDobj) from which the stream network
%           was derived.
%
% Output arguments
%
%     h     patch handle
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S);
%     plot3d(S,DEM)
%
% See also: STREAMobj, STREAMobj/plot, STREAMobj/plotdz
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 18. April, 2015


[x,y,z] = STREAMobj2XY(S,DEM);

baselevel = min(z);
baselevel = repmat(baselevel,size(z));
baselevel(isnan(z)) = nan;

n = numel(z);

x = [x;x];
y = [y;y];
z = [z;baselevel];

% create triangles
TIN = [bsxfun(@plus,(1:n)',[0 n n+1]);bsxfun(@plus,(1:n)',[0 n+1 1])];
TIN(any(TIN>size(x,1),2),:) = [];
I   = any(isnan(z(TIN)),2);
TIN(I,:) = [];

htemp = trisurf(TIN,x,y,z);
shading interp

if nargout == 1;
    h = htemp;
end
