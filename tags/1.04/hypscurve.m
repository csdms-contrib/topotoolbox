function varargout = hypscurve(dem,bins)

% plot hypsometric curve of a digital elevation model
%
% Syntax
%
%     hypscurve(dem)
%     hypscurve(dem,bins)
%     ax = hypscurve(...)
%     [rf,elev] = hypscurve(...)
%
% Description
%
%     A hypsometric curve is an empirical cumulative distribution function
%     of elevations in a catchment. hypscurv plots the hypsometric curve or
%     returns the relative frequencies of elevation values. Optionally, 
%     hypscurve bins the data in equally spaced containers. 
%
% Input
%
%     dem       elevation data
%     bins      number of bins
%
% Output
%
%     ax        axis handle
%     rf        relative frequencies (in %)
%     elev      elevation
%
% Example
%
%     load exampleDEM
%     hypscurve(dem,50)
%
% See also: hist
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 15. July, 2010

dem = dem(~isnan(dem));
dem = dem(:);

if nargin == 2;
    [n,elev] = hist(dem,bins); 
    n = flipud(n(:));
    elev = flipud(elev(:));
    n = cumsum(n);
    linestyle = 'o-';
else
    elev = sort(dem);
    elev = flipud(elev(:));
    
    % return unique entries only
    I = [true;diff(elev)~=0];
    elev = elev(I);
    n    = find(I);
    linestyle = '-';
 
end

% get relative frequencies
n = n./n(end) * 100;

% plot results
if nargout ~= 2;
    axis_handle = plot(n,elev,linestyle);
    axis xy
    xlabel('frequency [%]')
    ylabel('elevation [m]')
end

% prepare output
if nargout == 1
    varargout{1} = axis_handle;
elseif nargout == 2;
    varargout{1} = n;
    varargout{2} = elev;
end
    



