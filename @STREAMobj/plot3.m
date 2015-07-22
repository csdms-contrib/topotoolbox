function h = plot3(S,DEM,varargin)

% 3d-line plot of a STREAMobj
%
% Syntax
%
%     plot3(S,DEM)
%     plot3(S,z)
%     h = ...
%
% Description
%
%     The plot3 function displays a three-dimensional plot of a stream 
%     network S derived from a DEM.
%
% Input arguments
%
%     S      STREAMobj
%     DEM    digital elevation model (DEM)
%     z      node attribute list
%
% Output arguments
%
%     h      line handle
%
% See also: STREAMobj, STREAMobj/plot, STREAMobj/plot3d, STREAMobj/plotdz
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 25. June, 2015

narginchk(2,inf)

[x,y,z] = STREAMobj2XY(S,DEM);
htemp = plot3(x,y,z,varargin{:});

if nargout == 1
    h = htemp;
end
