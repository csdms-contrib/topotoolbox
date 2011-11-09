function demfs = fillsinks(dem,maxdepth)

% fill/remove pits, sinks or topographic depressions
%
% Syntax
%
%     demfs = fillsinks(dem)
%     demfs = fillsinks(dem,maxdepth)
%     demfs = fillsinks(dem,sinks)
%
% Description
%
%     fillsinks removes topographic depressions in a Digital Elevation 
%     Model (DEM). Use this function to enable a continuous flow towards 
%     the DEM edges. 
%
%     Sinks may, however, be closed basins or dolines and as such they are 
%     important features of DEMs. In order to account for such sinks, 
%     fillsinks allows you to specify a maximum depth of sinks, that 
%     will be filled, or to employ a logical grid (sinks) that is true
%     where sinks should remain (minima imposition).
%
% Input
%
%     dem       digital elevation model
%     maxdepth  positive scalar with maximum depth of sinks that will be
%               filled
%     sinks     logical matrix same size as dem with true elements
%               referring to sinks
%
% Output
%
%     demfs     digital elevation model with filled sinks
%
% Example
%
%     dem = peaks(100);
%     demfs = fillsinks(dem,5);
%     subplot(1,2,1);
%     surf(dem)
%     title('unfilled DEM')
%     subplot(1,2,2)
%     surf(demfs)
%     title('DEM with sinks with depth less than 5 filled')
%
%
%
% See also: IMFILL, IMRECONSTRUCT, IMIMPOSEMIN
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 6. July, 2011


% check input arguments
error(nargchk(1, 2, nargin))
if nargin == 1;
    
else
    if isscalar(maxdepth)
        if maxdepth<=0;
            error('TopoToolbox:incorrectinput',...
                    'maxdepth must not be zero or negative')
        end
    else
        if ~isequal(size(dem),size(maxdepth));
            error('TopoToolbox:incorrectinput',...
                    'sinks must have same size as the digital elevation model')
        end
        SINKS = maxdepth;
    end
        
end

% identify nans
Inan      = isnan(dem);
% set nans to -inf
dem(Inan) = -inf;

if nargin == 1;
    % fill depressions using imfill with an 8-neighborhood
    demfs      = imfill(dem,8,'holes');
    
elseif nargin==2 && isscalar(maxdepth);
    
    % create mask
    % complement image
    dem  = imcomplement(dem);
    
    % find regional maxima
    Imax = imregionalmax(dem,8);
    
    % create marker
    marker = -inf(size(dem));
    marker(Imax) = dem(Imax)-maxdepth;
    
    demfs  = imreconstruct(marker,double(dem));
    
    % are there any regional maxima that have been filled
    % to exactly maxdepth
    II  = demfs(Imax) <= dem(Imax)-maxdepth;
    
    
    % if yes, repeat morphological reconstruction with new marker
    if any(II)
        Imax = find(Imax);
        marker(Imax(II)) = dem(Imax(II));
        demfs  = imreconstruct(marker,double(dem));
    end
    
    % complement image again
    demfs = imcomplement(demfs);
else
    % minima imposition
    I = true(size(dem));
    I(2:end-1,2:end-1) = false;
    if any(Inan);
        I = (imdilate(Inan,ones(3)) & ~Inan) | I;
    end
    
    SINKS = SINKS | I;
    
    marker = -inf(size(dem));
    marker(SINKS) = -dem(SINKS);
    demfs = -imreconstruct(marker,-dem);
    
end

% nans in the dem are set to nan again
demfs(Inan) = nan;

