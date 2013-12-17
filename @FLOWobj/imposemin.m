function DEM = imposemin(FD,DEM,sl)

% minima imposition (carving) along drainage network
%
% Syntax
%
%     DEMc = imposemin(FD,DEM)
%     DEMc = imposemin(FD,DEM,sl)
%
% Description
%
%     imposemin carves a digital elevation model by downstream minima
%     imposition such that
%
%     DEM(ixc) = min(DEM(ix),DEM(ixc))
%
%     where ix is a linear index in DEM and ixc is the linear index of
%     the downward neighbor of ix. Depending on the breach length, long
%     flat channel section may be generated which can be avoided by
%     imposing a slight, minimum downward gradient sl such that
%
%     DEM(ixc) = min(DEM(ix)-dx*sl,DEM(ixc)) 
%
%     Take care to not choose a high gradient which will cause that
%     channels dip below the true land surface. 
%
% Input arguments
%
%     FD        flow direction object (FLOWObj)
%     DEM       digital elevation model (GRIDobj)
%     sl        minimum gradient [m/m] in downward direction (e.g. 0.001)
%
% Output arguments
%
%     DEMc      carved DEM
%
% 
%
% See also: FLOWobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 10. February, 2012


narginchk(2,3)
validatealignment(FD,DEM);
dem = DEM.Z;

if nargin == 2;    
    for r = 1:numel(FD.ix);
        dem(FD.ixc(r)) = min(dem(FD.ix(r)),dem(FD.ixc(r)));
    end
    
elseif nargin == 3;
    d = getdistance(FD.ix,FD.ixc,FD.size,FD.cellsize);
    d = d*sl;
    for r = 1:numel(FD.ix);
        dem(FD.ixc(r)) = min(dem(FD.ix(r))-d(r),dem(FD.ixc(r)));
    end
end

DEM.Z = dem;


