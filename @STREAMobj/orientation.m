function theta = orientation(S,varargin)

%ORIENTATION stream orientation
%
% Syntax
%
%     theta = orientation(S,pn,pv,...)
%
% Description
%
%     This function calculates the orientation of a stream for each node in
%     the network. Owing to the large variability of river directions, the
%     function smoothes the planform pattern of the network using
%     STREAMobj/smooth.
%
% Input arguments
%
%     S     STREAMobj
%     
%     Parameter name/value pairs pn,pv
%
%     'K'     smoothing factor (scalar > 0). Default = 100
%     'unit'  {'degrees'} or 'radians' or 'cart'. Degrees will be measured clockwise
%             from top, wheras radians will be measured anti-clockwise from
%             the x-axis. 'cart' returns the cartesian coordinates of the
%             directional vectors.
%
% Output arguments
%
%     theta   node-attribute list of angles of orientation 
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'mex',true,'preprocess','carve');
%     S  = STREAMobj(FD,'minarea',1000);
%     theta = orientation(S,'k',1000);
%     plotc(S,theta)
%     caxis([0 360])
%     colormap('hsv')
%
% See also: STREAMobj, SWATHobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 17. August, 2017


% Input checking and parsing
p = inputParser;
p.FunctionName = 'STREAMobj/orientation';
addParameter(p,'K',100,@(x) isscalar(x));
addParameter(p,'unit','deg');
parse(p,varargin{:});

unit = validatestring(p.Results.unit,{'radians','degrees','cart'},'STREAMobj/orientation');

y = smooth(S,S.y,'K',p.Results.K);
x = smooth(S,S.x,'K',p.Results.K);

dx = gradient(S,x);
dy = gradient(S,y);

switch unit
    case 'cart'
        theta = [dx dy];
        return
end

theta = cart2pol(dx,dy);

switch unit
    case 'degrees'
        theta = rad2deg(theta);
        theta = mod(-90-theta,360);
end