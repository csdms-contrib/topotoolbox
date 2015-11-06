function h = plotc(S,DEM)

% plot a colored stream network
%
% Syntax
%      
%     plotc(S,DEM)
%     plotc(S,nal)
%     h = plotc(...)
%
% Description
%
%     plotc plots the planform stream network with additional coloring
%     obtained from a GRIDobj or a node attribute list.
%
% Input arguments
%
%     S      STREAMobj
%     DEM    GRIDobj
%     nal    node attribute list
%
% Output arguments
%
%     h      handle to surface object
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','carve');
%     S   = STREAMobj(FD,'minarea',1000);
%     subplot(1,2,1);
%     imageschs(DEM,DEM,'colormap','gray','colorbar',false)
%     hold on
%     plotc(S,DEM)
%     colorbar
%     subplot(1,2,2);
%     imageschs(DEM,DEM,'colormap','gray','colorbar',false)
%     hold on
%     s = streamorder(S);
%     plotc(S,s)
%     colorbar     
%
%
% See also: STREAMobj/plot, STREAMobj/plot3, STREAMobj/plot3d,
%           STREAMobj/plotdz
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 6. November, 2015


[x,y,c] = STREAMobj2XY(S,DEM);
z = x*0;
ht = surface([x x],[y y],[z z],[c c],...
        'facecolor','none',...
        'edgecolor','flat',...
        'linewidth',2);
if nargout == 1;
    h = ht;
end
