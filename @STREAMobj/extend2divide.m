function S = extend2divide(S,FD)

%EXTENT2DIVIDE Grow STREAMobj upstream to extend to the divide
%
% Syntax
%
%     S2 = extend2divide(S,FD)
%
% Description
%
%     extend2divide takes a stream network (STREAMobj) S and a FLOWobj FD
%     and creates a new network S2 which comprises S and single channels 
%     expanding from each channelhead in S to the divide. Since there are 
%     multiple paths available, the function chooses the one along the
%     highest values of flow accumulation.
%
% Input arguments
%
%     S      STREAMobj
%     FD     FLOWobj (S must have been derived based on FD)
%
% Output arguments
%
%     S2     STREAMobj
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD   = FLOWobj(DEM);
%     S = STREAMobj(FD,'minarea',1000);
%     S2 = extend2divide(S,FD);
%     plot(S2)
%     hold on
%     plot(S)
%
%  
% See also: STREAMobj/modify, STREAMobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 5. August, 2021

% create a copy of FD
FDcopy = FD;
% find channelheads
ix     = streampoi(S,'ch','ix');
% Determine contributing areas of channelheads
D  = dependencemap(FD,ix);
% Clip the FLOWobj to dependence map
I  = D.Z(FD.ix) & D.Z(FD.ixc);
FD.ix = FD.ix(I);
FD.ixc = FD.ixc(I);
% Calculate flow accumulation
A  = flowacc(FD);
% Flip flow directions
FD = flipdir(FD);
% Recalculate flow accumulation, this time in reverse direction and with 
% previously calculated accumulation as weight,
A  = flowacc(FD,A);
% Use these values as weights for the multiple (reversed) flow direction
FD.fraction = A.Z(FD.ixc);
% Calculate single from multiple flow directions
FD = multi2single(FD);
% Use channelheads as seeds
W  = GRIDobj(A);
W.Z(ix) = 1;
% Where do these seeds end up
A  = flowacc(FD,W);
[~,ix] = drainagebasins(FD);
I = A.Z(ix) == 1;
ix = ix(I);
FD = FDcopy;
S2  = STREAMobj(FD,'ch',ix);
S  = modify(S2,'upstreamto',streampoi(S,'outl','ix'));


