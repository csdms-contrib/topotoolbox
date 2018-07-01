function DEM = imposemin(S,DEM,sl)

%IMPOSEMIN minima imposition (carving) along stream network
%
% Syntax
%
%     DEMc = imposemin(S,DEM)
%     DEMc = imposemin(S,DEM,sl)
%     zc   = imposemin(S,z)
%     zc   = imposemin(S,z,sl)
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
%     z         node attribute list of elevation values 
%     sl        minimum gradient [m/m] in downward direction (e.g. 0.001)
%
% Output arguments
%
%     DEMc      carved DEM
%     zc        node attribute list of elevation values
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','c');
%     S  = STREAMobj(FD,A>1000);
%     S  = klargestconncomps(trunk(S));
%     plotdz(S,DEM)
%     hold on
%     DEMc = imposemin(S,DEM);
%     plotdz(S,DEMc)
%
% See also: FLOWobj/imposemin
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 17. August, 2017


narginchk(2,3)

if isa(DEM,'GRIDobj')
    inp = 'GRIDobj';
    validatealignment(S,DEM);
    z = getnal(S,DEM);
elseif isnal(S,DEM);
    inp = 'nal';
    z = DEM;
else
    error('Imcompatible format of second input argument')
end
    
if nargin == 2;
    for r = 1:numel(S.ix);
        z(S.ixc(r)) = min(z(S.ix(r)),z(S.ixc(r)));
    end    
elseif nargin == 3;
    
    d = sqrt((S.x(S.ix)-S.x(S.ixc)).^2 + (S.y(S.ix)-S.y(S.ixc)).^2); 
    d = d*sl;
    for r = 1:numel(S.ix);
        z(S.ixc(r)) = min(z(S.ix(r))-d(r),z(S.ixc(r)));
    end
end

switch inp
    case 'GRIDobj'
        DEM.Z(S.IXgrid) = z;
    case 'nal'
        DEM = z;
end


