function R = roughness(dem,type,ks)

% terrain ruggedness, position and roughness indices of DEMs
%
% Syntax
%
%     R = roughness(dem,type,ks)
%
% Description
%
%     roughness calculates various metrics for locally quantifying terrain
%     roughness. Currently the function supports three different roughness
%     indices. See description below.
%
% Input
%
%     dem   digital elevation model
%     type  roughness metrics type 
%           'tri' topographic ruggedness index (Riley et. al. 1999)
%                 (corresponds to 2d standard deviation filter
%           'tpi' topographic position index
%                 (the difference between elevation in a pixel and mean
%                 elevation of its surrounding pixels)
%           'roughness' Roughness is the the largest inter-cell difference 
%                 of a central pixel and its surrounding cell (=default)
%     ks    kernel size (default = [3 3])
% 
% Output
%     
%     R     roughness index
%
% References
% 
%     Riley, S. J., DeGloria, S. D., Elliot, R. 1999:  A terrain ruggedness 
%     index that quantifies topographic heterogeneity. Intermountain
%     Journal of Sciences, 5, 23-27.
%
% See also: stdfilt
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 24. February, 2010

if nargin == 1;
    type = 'roughness';
    ks   = [3 3];
elseif nargin == 2;
    ks   = [3 3];
end

% check if kernel has right size
if numel(ks) ~= 2;
    error('TopoToolbox:incorrectinput',...
        'ks must be a two element vector');
end

% check if kernel rows and cols are uneven
if any(mod(ks,2),1) ~= 1
    error('TopoToolbox:incorrectinput',...
        'the kernel must have uneven number of rows and cols (e.g. [3 3] or [5 5])');
end

% pad for some roughness indices
switch lower(type)
    case {'tpi'}       
        
        padrc = floor(ks/2);
        dem = padarray(dem,padrc,'symmetric');
        depad = true;
    otherwise
        depad = false;
          
end

I   = isnan(dem);
[L,L] = bwdist(~I);
dem     = dem(L);

switch lower(type)
    case 'tri'
        % terrain ruggedness index
        % root mean squared difference between center pixel and surrounding
        % pixels
        % excerpt from ArcGIS-Help
        % The topographic ruggedness index (TRI) is a measurement developed
        % by Riley, et al. (1999) to express the amount of elevation 
        % difference between adjacent cells of a digital elevation grid. 
        % The process essentially calculates the difference in elevation 
        % values from a center cell and the eight cells immediately 
        % surrounding it. Then it squares each of the eight elevation 
        % difference values to make them all positive and averages the 
        % squares. The topographic ruggedness index is then derived by 
        % taking the square root of this average, and corresponds to 
        % average elevation change between any point on a grid and it’s 
        % surrounding area.
        %
        % I don't have access to the article, but it seems to me that this
        % index is nothing but a moving standard deviation
        
        kernel = ones(ks);
        R      = stdfilt(dem,kernel);
        
    case 'tpi'
        % Topographic Position Index
        % difference between a central pixel and the mean of its
        % surrounding cells
        kernel = ones(ks);
        kernel(ceil(ks(1)/2),ceil(ks(2)/2)) = 0;
        % normalize kernel
        kernel = kernel/sum(kernel(:));
        
        % convolute
        R = conv2(dem,kernel,'valid');
        R = dem(padrc(1)+1:end-padrc(1),padrc(2)+1:end-padrc(2))-R;      
        
    case 'roughness'
        % Roughness
        % Roughness is the the largest inter-cell difference of a central
        % pixel and its surrounding cell
        kernel = ones(3);        
        R = max(imdilate(dem,kernel)-dem,dem-imerode(dem,kernel));
        
        if nargin == 3
        warning('TopoToolbox:incorrectinput',...
            'the kernel size does not apply to the roughness index.')
        end
        
end
        
if depad
    R(I(padrc(1)+1:end-padrc(1),padrc(2)+1:end-padrc(2))) = nan;
else
    R(I) = nan;
end