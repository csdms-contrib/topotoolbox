function [FD,S,DA,dx_ik] =updateDrainDir_hylands(FlowDir,H1,X,Y,boundaryNodesAll,nb,FD)
% Function to update the drainage network
%
% Syntax
%
%       [FD,C,S,i,k,p,DA,A,dx_ik] = updateDrainDir(H1,BORDER,W,p,X,Y,R)
%
% Description
%
%       Function to update the drainage network 
%
% Input
%
%       H1        DEM (digital elevation model) (GRIDobj)
%       BORDER    (GRIDobj) produced with getBORDER function
%       W         Weights
%       p         structure array with parameter definitions (see ttlemset)
%       X         gridded X distance Y         gridded Y distance R
%       Randomized matrix with dimensions of DEM.Z
%
%
% Output
%
%       FD        Flow Direction, FLOWobj
%       C         River cells
%       S         STREAMobj
%       i         Giver
%       k         Receiver
%       p         updated structure array with parameter definitions (see ttlemset)
%       A         Velocity field of advective stream power law
%       dx_ik     distance between giver and receiver
%       kk        Receiver of receiver
%       ii        Giver of giver
%       dx_centered distance between giver of giver and receiver
%
% Example
%
%
% See also: HYLANDS, HYLANDS_set
%
% Authors:  Benjamin Campforts (benjamin.campforts@kuleuven.be)
%           Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
%
% Date:     15. March, 2020
%
%------------------------------ References -------------------------------%
% HyLands: Campforts B., Shobe M.C., et al. : HyLands 1.0: a Hybrid
% Landscape evolution model to simulate the impact of landslides and
% landslide-derived sediment on landscape evolution. Discussion paper in
% Geoscientific Model Development,
% https://geoscientific-model-development.net
%
% TTLEM: Campforts B., Schwanghart W., Govers G. (2016): TTLEM 1.0 : a
% numerical package for accurate simulationsng  of transient landscape
% evolution in MATLAB. Discussion paper in GMD.
%
%-------------------------------------------------------------------------%


switch FlowDir
    case 'single'
        if isempty(FD)
            FD = FLOWobj(H1,'mex',true);
        end
    case 'multi'
        FD_mult = multiFLOWobj(H1,'multi',nb);
        FD=FLOWobj;
        FD.size=FD_mult.size;
        FD.type=FD_mult.type;
        FD.ix=FD_mult.ix;
        FD.ixc=FD_mult.ixc;
        FD.fraction=FD_mult.fraction;
        FD.cellsize=FD_mult.cellsize;
        FD.refmat=FD_mult.refmat;
        clear FD_mult
    case 'DInf'
        FD = FLOWobj(H1,'DInf');
end


dx2=FD.cellsize.^2;
W0=GRIDobj(H1)+1;
W0.Z(boundaryNodesAll)=0;
DA  = flowacc(FD,W0); 
A = GRIDobj(H1);
A.Z = flowacc_mex(FD.ix,FD.ixc,FD.size);
DA=DA.*dx2;
S = [];
switch FlowDir
    case {'DInf','multi'}
    otherwise
        FD.fraction=ones(size(FD.ix));
end
dx_ik = double(sqrt((X(FD.ix)-X(FD.ixc)).^2 + (Y(FD.ix)-Y(FD.ixc)).^2));
