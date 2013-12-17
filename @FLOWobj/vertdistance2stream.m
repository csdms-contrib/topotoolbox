function DZ = vertdistance2stream(FD,S,DEM)

% vertical distance to streams
%
% Syntax
%
%     DZ = vertdistance2stream(FD,S,DEM)
%
% Description
%
%     vertdistance2stream calculates the vertical distance of each cell in
%     a digital elevation model DEM to the nearest stream cell in S along
%     the flow path in FD.
%
% Input arguments
%
%     FD    instance of FLOWobj
%     DEM   digital elevation model (class: GRIDobj)
%     S     stream network (class: STREAMobj)
%
% Output arguments
%
%     DZ    vertical distance to streams (class: GRIDobj)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','c');
%     S = STREAMobj(FD,'minarea',1e6,'unit','m');
%     DZ = vertdistance2stream(FD,S,DEM);
%     imageschs(DEM,DZ)
%     hold on
%     plot(S,'k','LineWidth',2)
% 
%
% See also: FLOWobj, FLOWobj/flowdistance, GRIDobj, STREAMobj
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 25. January, 2013

narginchk(3,3)

validatealignment(S,DEM);
validatealignment(FD,DEM);

DZ = DEM;
DZ.Z = -inf(DEM.size);
DZ.Z(S.IXgrid) = DEM.Z(S.IXgrid);

for r = numel(FD.ix):-1:1
    DZ.Z(FD.ix(r)) = max(DZ.Z(FD.ix(r)),DZ.Z(FD.ixc(r)));
end
DZ = DEM-DZ;
