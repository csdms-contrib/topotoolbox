function S = split(S,varargin)

% split drainage network at predefined locations
%
% Syntax
%
%     S = split(S)
%     S = split(S,ix)
%
% Description
%
%     split disconnects the channel network at predefined locations. With
%     only the STREAMobj as input argument, the stream network will be
%     split at all confluences. Otherwise, the network will be split at ix 
%     which is a vector of linear indices into the DEM from which S was 
%     derived. 
%
% Input arguments
%
%     S     STREAMobj
%     ix    linear index into DEM. The indexed pixels must be located on 
%           the channel network, e.g. ismember(ix,S.IXgrid).
% 
% Output arguments
%
%     Ss    STREAMobj
%
% Example 1: Split at confluences
%     
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1e6,'unit','map');
%     S2 = split(S);
%     plot(S)
%     hold on
%     plot(S2,'k')
%     % zoom in to see that the stream network has been split
%
% Example 2: Split at random locations
%
%     S2 = split(S,randlocs(S,100));
%     plotc(S2,S2.distance)
%
% See also: STREAMobj, STREAMobj/modify, STREAMobj/randlocs, 
%           STREAMobj/STREAMobj2cell
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 19. May, 2016


narginchk(1,2)
if nargin == 1;
    V = streampoi(S,'confluences','logical');
    I = V(S.ixc);
else
    I = ismember(S.IXgrid,varargin{1});
    I = I(S.ixc);
end

S.ix(I) = [];
S.ixc(I) = [];