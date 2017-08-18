function [lat,lon,varargout] = STREAMobj2latlon(S,varargin)

%STREAMOBJ2LATLON convert instance of STREAMobj to NaN-separated geographic coordinates
%
% Syntax
%
%     [lat,lon] = STREAMobj2latlon(S)
%     [lat,lon,a,b,...] = STREAMobj2latlon(S,A,B,...)
%
% Description
%
%     STREAMobj2latlon returns the NaN-punctuated latitude and longitude
%     vectors which can be used for easy plotting using the plot function.
%     With additional input arguments that are instances of GRIDobj,
%     additional vectors are generated with the respective values of the
%     grids A, B, etc. at the node locations of the stream network S.
%
%     Note that this function *requires the mapping toolbox*.
%
% Input arguments
%
%     S       streams (class STREAMobj)
%     A,B,... grids (class GRIDobj) or node attributes (e.g. as returned by
%             the function STREAMobj/streamorder
%
% Output arguments
%
%     lat,lon  coordinate vectors 
%     a,b,...  grid values at locations specified by x and y.
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     [lat,lon] = STREAMobj2latlon(S);
%     plot(lon,lat)
%
% See also: STREAMobj/STREAMobj2XY
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 20. February, 2015

nr = numel(varargin);
c  = cell(1,nr);
[x,y,c{:}] = STREAMobj2XY(S,varargin{:});
if isempty(S.georef)
    error('The STREAMobj contains no georeferencing information');
end
[lat,lon] = minvtran(S.georef.mstruct,x,y);
varargout = c;