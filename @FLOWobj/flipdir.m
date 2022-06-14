function FD = flipdir(FD,varargin)

%FLIPDIR Flip direction of flow
%
% Syntax
%
%     FDm = flipdir(FD)
%
% Description
%
%     flipdir changes the direction of flow in a single or multi FLOWobj so
%     that flow moves upstream. In upstream direction, single FLOWobjs
%     become divergent and thus flipdir always returns a FLOWobj with
%     multiple flow directions.
%
% Input arguments
%
%     FD     FLOWobj
%
% Output arguments
%
%     FDm    flipped FLOWobj (multi)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(fillsinks(DEM),'multi');
%     [FD,S]  = multi2single(FD,'minarea',500);
%     FDf = flipdir(FD);
%     H = ~STREAMobj2GRIDobj(S);
%     A = flowacc(FDf,H);
%     imageschs(DEM,min(A,1000))
%
%
% See also: FLOWobj, FLOWobj/multi2single
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 24. December, 2021

A = flowacc(FD);

FD.type = 'multi';
ix = FD.ix;
FD.ix = FD.ixc;
FD.ixc = ix;

% flip ordering
FD.ix  = FD.ix(end:-1:1);
FD.ixc = FD.ixc(end:-1:1);
if isempty(FD.fraction)
    FD.fraction = ones(size(FD.ix));
else
    
    FD.fraction = A.Z(FD.ixc);
%     FD.fraction = FD.fraction(end:-1:1);
end
FD = multi_normalize(FD);


