function DZ = vertdistance2stream(FD,S,DEM)

%VERTDISTANCE2STREAM vertical distance to streams
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
% See also: FLOWobj, FLOWobj/flowdistance, FLOWobj/mapfromnalGRIDobj, 
%           STREAMobj
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 4. March, 2016


% 4/3/2016: the function now makes copies of FD.ix and FD.ixc (see 
% FLOWobj/flowacc

narginchk(3,3)

validatealignment(S,DEM);
validatealignment(FD,DEM);

DZ = DEM;
DZ.Z = -inf(DEM.size);
DZ.Z(S.IXgrid) = DEM.Z(S.IXgrid);

ix = FD.ix;
ixc = FD.ixc;
for r = numel(ix):-1:1
    DZ.Z(ix(r)) = max(DZ.Z(ix(r)),DZ.Z(ixc(r)));
end
DZ = DEM-DZ;
