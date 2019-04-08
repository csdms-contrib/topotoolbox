function zi = interp(S,d,z,varargin)

%INTERP interpolate data on STREAMobj (single river only)
%
% Syntax
%
%     zi = interp(S,d,z)
%     zi = interp(S,d,z,'method','extrapolation')
%
% Description
%
%     
% Input arguments
%
%     S      STREAMobj with one channelhead
%     d      distance values
%     z      variable values (same size as d)
%     
%     'method'  see interp1 for details
%     'extrapolation'  see interp1 for details
%
% Output arguments
%
%     zi     node-attribute list with interpolated values
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','c');
%     S   = STREAMobj(FD,'minarea',1000);
%     S   = trunk(klargestconncomps(S));
%     zi = interp(S,[10000 12000 20000 39000 40000],[5 10 20 6 10],...
%                 'spline',10);
%     stackedplotdz(S,{DEM zi})
%     
% See also: STREAMobj, STREAMobj/trunk, demo_modifystreamnet
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 18. February, 2019


p = inputParser;
addRequired(p,'S',@(x) numel(streampoi(x,'channelheads','ix')) == 1);
addRequired(p,'d')
addRequired(p,'z',@(x) isequal(size(x),size(d)));
addOptional(p,'method','linear',@(x) ischar(x) || isstring(x));
addOptional(p,'extrapolation','extrap',@(x) ischar(x) || isstring(x) || isscalar(x));
parse(p,S,d,z,varargin{:});
S   = p.Results.S;

[d,ix] = sort(d,'ascend');
z = z(ix);

zi = interp1(d,z,S.distance,p.Results.method,p.Results.extrapolation);
