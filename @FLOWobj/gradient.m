function G = gradient(FD,DEM)

%GRADIENT gradient along flow direction
%
% Syntax 
%
%     G = gradient(FD,DEM)
%
% Description
%
%     gradient (slope) along flow direction not necessarily follows the
%     steepest descend when carving or filling were chosen as preprocessing
%     when deriving the FLOWobj. FLOWobj2gradient returns an instance of
%     GRIDobj that contains the gradient along the flow paths.
%
%     If FLOWobj has been derived using a multiple or Dinf flow direction
%     algorithm, then G is calculated as the weighted mean gradient.
%
% Input arguments
%
%     FD     FLOWobj
%     DEM    digital elevation model (GRIDobj)
%
% Output arguments
%
%     G      along-flow gradient
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','c');
%     G  = gradient(FD,DEM);
%     imageschs(DEM,G)
%
% See also: GRIDobj/gradient8, GRIDobj/arcslope, STREAMobj/gradient
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 27. March, 2019

validatealignment(FD,DEM)
d = getdistance(FD.ix,FD.ixc,FD.size,FD.cellsize);
G = DEM;
switch FD.type
    case 'single'
        G.Z = zeros(FD.size,class(DEM.Z));
        G.Z(FD.ix) = (DEM.Z(FD.ix)-DEM.Z(FD.ixc))./d;
    otherwise
        d   = FD.fraction .* (DEM.Z(FD.ix) - DEM.Z(FD.ixc))./d;
        g   = accumarray(FD.ix,d,[prod(DEM.size) 1],@sum,zeros(1,'like',DEM.Z),false);
        G.Z = reshape(g,DEM.size);
end
        