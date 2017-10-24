function S = removeedgeeffects(S,FD,DEM)

%REMOVEEDGEEFFECTS remove potential edge effects
%
% Syntax
%
%     S2 = removeedgeeffects(S,FD)
%     S2 = removeedgeeffects(S,FD,DEM)
%
% Description
%
%     REMOVEEDGEEFFECTS removes nodes in the stream network S that have 
%     boundary pixels in their catchment area. Removing these nodes ensures
%     that the stream network contains only nodes whose upstream area is
%     completely covered by the DEM.
%
%     removeedgeeffects(S,FD) removes only those pixels that are downstream
%     of pixels that are located on the edges of the grid GRIDobj(FD).
%
%     removeedgeeffects(S,FD,DEM) accounts for nan-pixels in the DEM.
%
% Input arguments
%
%     S     STREAMobj
%     FD    FLOWobj
%     DEM   Digital elevation model (GRIDobj). 
%
% Output arguments
%
%     S2    trimmed STREAMobj
%
% 
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM);
%     S  = STREAMobj(FD,'minarea',1e6,'unit','map');
%     plot(S)
%     S  = removeedgeeffects(S,FD);
%     hold on
%     plot(S)
%
% See also: STREAMobj/mnoptim, STREAMobj/chitransform
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 24. October, 2017

narginchk(2,3)

if nargin == 2;
    I  = GRIDobj(FD,'logical');
    I.Z(:,:) = true;
    I.Z(2:end-1,2:end-1) = false;
else

    II = isnan(DEM);
    nans = any(II);
    I  = GRIDobj(DEM,'logical');
    I.Z(:,:) = true;
    I.Z(2:end-1,2:end-1) = false;
    
    if nans
        II = dilate(II,ones(3));
        I  = I | II;
        clear II
    end
end


I  = influencemap(FD,I);
S = modify(S,'upstreamto',~I);