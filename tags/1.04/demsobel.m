function dems = demsobel(dem,type)

% edge detection using sobel filter 
%
% Syntax
%
%     dems = demsobel(dem)
%     dems = demsobel(dem,type)
%
% Description
% 
%     demsobel
%
% Input
%
%     dem       digital elevation model
%     type      'sobel' or 'scharr'
%
% Output
%
%     dems      smoothed digital elevation model
%
% Example
%     
%     % highlight ridges and valleys by substracting
%     % the original digital elevation model from
%     % the smoothed one.
%
%     load exampleDEM
%     dems = demsmooth(dem,11);
%     surf(X,Y,dem,dems-dem); 
%     title('highlight ridges and valleys')
% 
% See also: CONV2, FILTER2, BWDIST
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 15. March, 2009


% ..............................................
% check input arguments
if nargin == 1;
    type = 'sobel';
    
elseif nargin == 2;
else
    error('wrong number of inputs')
end



% ..............................................
% are there any nans?
inan    = isnan(dem);
flagnan = any(inan);

if flagnan;
    % if there are nans pad array with nans and use the bwdist function to
    % do a nearest neighbor interpolation and extrapolation
    dem     = padarray(dem,padsize,nan);
    inanpad = isnan(dem);
    [L,L]   = bwdist(~inanpad);
    dem     = dem(L);
else
    % if no nans use the usual padarray function with the replicate option
    dem     = padarray(dem,[1 1],'replicate');
end


switch lower(type)
    case 'sobel'
        ky = [1 2 1; 0 0 0; -1 -2 -1];
        kx = ky';
    case 'scharr'
        ky = [3 10 3; 0 0 0; -3 -10 -3];
        kx = ky';
    otherwise
        error('unknown operator')
end

dems = hypot(conv2(dem,ky,'valid'),conv2(dem,kx,'valid'));

if flagnan;
    dems(inan) = nan;
end



