function preprocesstool(DEM)

% interactive tool for hydrological conditioning
%
% Syntax
%
%     preprocesstool(DEM)
%
% Description
%
%     preprocesstool starts an interactive tool that enables the
%     hydrological conditioning of a DEM using sink filling and carving. 
%
% Input arguments
%
%     DEM    digital elevation model (class: GRIDobj)
%
% See also: FLOWobj/imposemin, GRIDobj/fillsinks
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 17. April, 2013

narginchk(1,1)
validateattributes(DEM,{'GRIDobj'},{},'preprocesstool','DEM',1)
gui_preprocess(DEM);