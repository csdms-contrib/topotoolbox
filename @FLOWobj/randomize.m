function FD = randomize(FD)

%RANDOMIZE randomize multiple flow directions
%
% Syntax
%
%     FDr = randomize(FDm)
%
% Description
%
%     randomize takes an instance of a multiple flowdirections FLOWobj and
%     randomizes the flow directions. This removes the slope dependency of
%     the fraction that each cell transfers to its downstream neighbors. 
%
% Input arguments
%
%     FDm     FLOWobj (multiple flow directions)
%
% Output arguments
%
%     FDr     FLOWobj (randomized multiple flow directions
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     DEM = filter(DEM,'mean',[21 21]);
%     DEM = fillsinks(DEM);
% 
%     ax(1) = subplot(1,2,1);
%     FD  = FLOWobj(DEM,'multi');
%     A = flowacc(FD);
%     imageschs(DEM,log(A),'colormap',flowcolor)
% 
%     ax(2) = subplot(1,2,2);
%     FDr = randomize(FD);
%     Ar = flowacc(FDr);
%     imageschs(DEM,log(Ar),'colormap',flowcolor)
% 
%     linkaxes(ax,'xy')
%     % Now zoom in to explore differences in flow patterns
%
%
% See also: FLOWobj, FLOWobj/multi2single, FLOWobj/flowconvergence
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 23. May, 2019

% If FD has single flow directions, do nothing
if ~ismulti(FD)
    return
end

% randomize the flow fraction using a cont. uniform distribution
FD.fraction = rand(size(FD.fraction));

% normalize so that the sum in each fraction is one
FD = multi_normalize(FD);

% Normalize code
% [~,~,ix] = unique(FD.ix);
% s = accumarray(ix,FD.fraction);
% s = s(ix);
% FD.fraction = FD.fraction./s;
