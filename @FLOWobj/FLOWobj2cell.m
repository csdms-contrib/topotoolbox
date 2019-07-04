function [CFD,D,A,cfix] = FLOWobj2cell(FD,IX)

%FLOWOBJ2CELL return cell array of FLOWobjs for individual drainage basins
%
% Syntax
%
%     CF = FLOWobj2cell(FD)
%     [CF,D,a] = FLOWobj2cell(FD)
%     [CF,D,a,cfix] = FLOWobj2cell(FD,IX)
%
% Description
%
%     FLOWobj2cell derives a cell array of FLOWobjs for each individual
%     drainage basin. This may be required if computations on individual
%     drainage basins should be run in parallel.
%
% Input arguments
%
%     FD    FLOWobj
%     IX    linear index in GRIDobj from which FD was derived. If IX is
%           provided, the function returns a forth output argument that
%           indicates in which of the flow networks (FLOWobjs) in the
%           output cell array IX is located.
%
% Output arguments
%
%     CF    cell array of FLOWobjs
%     D     GRIDobj with drainage basins
%     a     area of drainage basins (unit = nr of pixels)
%     cfix  (see input argument IX)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','c');
%     [CFD,~,a] = FLOWobj2cell(FD);
%     [~,ix] = max(a);
%     A = flowacc(CFD{ix});
%     imageschs(DEM,log(A))
%
% See also: FLOWobj, STREAMobj/STREAMobj2cell
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 9. January, 2015

D   = drainagebasins(FD);
nrd = max(D);
dix = D.Z(FD.ix);
[dix,I] = sort(dix,'ascend');

if nargin == 2
    cfix = D.Z(IX);
end

if nargin == 1 && nargout == 4
    error('a forth output argument is only allowed for two input arguments');
end



FDcopy = FD;
FDcopy.ix = [];
FDcopy.ixc = [];

IX = (1:numel(FD.ix))';

CFD   = accumarray(dix,IX(I),[nrd 1],@distribute2FLOWobj);
if nargout > 2
    A     = cellfun(@(FD) numel(FD.ix)+1,CFD);
end

function FO = distribute2FLOWobj(IX)

FO     = FDcopy;
FO.ix  = FD.ix(IX);
FO.ixc = FD.ixc(IX);
FO     = {FO};
end
end






