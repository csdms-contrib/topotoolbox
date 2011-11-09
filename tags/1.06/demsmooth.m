function dems = demsmooth(dem,ws,varargin)

% moving average filter for digital elevation models 
%
% Syntax
%
%     dems = demsmooth(dem,ws)
%
% Description
% 
%     demsmooth is a mean filter with window size ws. ws can be an integer 
%     scalar that defines the edge length of a square kernel or a two 
%     element vector defining the size of the kernel (e.g. ws = [3 5]).
%
%     demsmooth handles NaNs using a nearest neighbor interpolation
%
% Input
%
%     dem       digital elevation model
%     ws        kernel size (default: [3 3])
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
% Date: 6. July, 2011


% ..............................................
% check input arguments
if nargin == 1;
    ws = [3 3];
elseif nargin > 2;
    error('wrong number of inputs')
end



if isscalar(ws);
    ws = [ws ws];
end


% ..............................................
% are there any nans?
inan    = isnan(dem);
flagnan = any(inan);

padsize = ceil(ws/2);
if flagnan;
    % if there are nans pad array with nans and use the bwdist function to
    % do a nearest neighbor interpolation and extrapolation
    dem     = padarray(dem,padsize,nan);
    inanpad = isnan(dem);
    [L,L]   = bwdist(~inanpad); %#ok<ASGLU>
    dem     = dem(L);
else
    % if no nans use the usual padarray function with the replicate option
    dem     = padarray(dem,padsize,'replicate');
end

% ..............................................
% create the filter matrix

% usual mean filter
W = ones(ws);

% IDW filter (see Neson and Jones 1995)
% if any(mod(ws,2) == 0)    
%     W = false(ws);
%     W(ceil(ws(1)/2),ceil(ws(1)/2)) = true;
%     W = bwdist(W);
%     W(ceil(ws(1)/2),ceil(ws(1)/2)) = 1;
% end


dems = conv2(dem,W./sum(W(:)),'same');

% ..............................................
% depad array
dems = dems(padsize(1)+1:end-padsize(1),padsize(2)+1:end-padsize(2));

% and set nans at their previous position
if flagnan
    dems(inan) = nan;
end



