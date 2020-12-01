function h = histogram(varargin)

%HISTOGRAM Plot frequency distribution of values in GRIDobj
%
% Syntax
%
%     histogram(DEM)
%     histogram(DEM,...)
%     h = ...
%
% Description
%
%     HISTOGRAM overloads MATLAB's histogram function for GRIDobjs. Please
%     read the help for this function for details.
%
% Examples
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     histogram(DEM)
%
%  
% See also: GRIDobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 1. December, 2020

ixarg = cellfun(@(x) isa(x,'GRIDobj'),varargin,'UniformOutput',true);
DEM = varargin{ixarg};
z   = DEM.Z;
if isfloat(z)
    z(isnan(z) | isinf(z)) = [];
end

z   = z(:);
varargin{ixarg} = z;

htemp = histogram(varargin{:});
if nargout == 1
    h = htemp;
end

