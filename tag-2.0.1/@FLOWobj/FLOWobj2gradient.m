function G = FLOWobj2gradient(FD,DEM)

% gradient along flow direction
%
% Syntax 
%
%     G = FLOWobj2gradient(FD,DEM)
%
% Description
%
%     gradient (slope) along flow direction not necessarily follows the
%     steepest descend when carving or filling were chosen as preprocessing
%     when deriving the FLOWobj. FLOWobj2gradient returns an instance of
%     GRIDobj that contains the gradient along the flow paths.
%
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 23. February, 2013

validatealignment(FD,DEM)
d = getdistance(FD.ix,FD.ixc,FD.size,FD.cellsize);
G = DEM;
G.Z = zeros(FD.size,class(DEM.Z));
G.Z(FD.ix) = (DEM.Z(FD.ix)-DEM.Z(FD.ixc))./d;