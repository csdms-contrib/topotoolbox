function DEM = imposemin(S,DEM,sl)

% minima imposition (carving) along stream network
%
% Syntax
%
%     DEMc = imposemin(S,DEM)
%     DEMc = imposemin(S,DEM,sl)
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
%     STREAMobj/imposemin works like FLOWobj/imposemin except that minima
%     imposition only affects the stream network.
%
% Input arguments
%
%     S         stream object (FLOWObj)
%     DEM       digital elevation model (GRIDobj)
%     sl        minimum gradient [m/m] in downward direction (e.g. 0.001)
%
% Output arguments
%
%     DEMc      carved DEM
%
% 
%
% See also: FLOWobj/imposemin
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 6. June, 2014


narginchk(2,3)
validatealignment(S,DEM);
dem = DEM.Z;

if nargin == 2;    
    for r = 1:numel(S.ix);
        dem(S.IXgrid(S.ixc(r))) = min(dem(S.IXgrid(S.ix(r))),dem(S.IXgrid(S.ixc(r))));
    end
    
elseif nargin == 3;
    d = sqrt((S.x(S.ix)-S.x(S.ixc)).^2 + (S.y(S.ix)-S.y(S.ixc)).^2); 
    d = d*sl;
    for r = 1:numel(FD.ix);
        dem(S.IXgrid(S.ixc(r))) = min(dem(S.IXgrid(S.ix(r)))-d(r),dem(S.IXgrid(S.ixc(r))));
    end
end

DEM.Z = dem;


