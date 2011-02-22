function [flats,Imr] = identifyflats(dem)

% identify flat terrain in a digital elevation model
%
% Syntax
%
%     [I,SILLS] = identifyflats(dem)
%
% Description
%
%     identifyflats returns a logical matrix that is true for cells
%     indicating flat terrain. flat terrain cells are defined as cells that
%     do not have a downward neighboring cell. The second output argument 
%     contains a logical matrix that is true for sill cells. Sill cells are
%     pixels in the DEM where flat regions spill over into lower terrain.   
%
% Input
%
%     dem        digital elevation model
%    
% Output
% 
%     I          logical matrix (true cells indicate flat terrain)
%
% Example
%
%     dem = peaks(300);
%     dem = fillsinks(dem);
%     I = identifyflats(dem);
%     surf(dem,+I)
%     shading flat; camlight
%
% 
% See also: ROUTEFLATS, CROSSFLATS
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 6. September, 2010



% handle NaNs
log_nans = isnan(dem);
if any(log_nans(:));
    flag_nans = true;
    dem(log_nans) = -inf;
else
    flag_nans = false;
end
nhood = ones(3);


% identify flats
% flats: logical matrix with true where cells don't have lower neighbors
if flag_nans
    flats = imerode(dem,nhood) == dem & ~log_nans;
else
    flats = imerode(dem,nhood) == dem;
end

% remove flats at the border
flats(1:end,[1 end])  = false;
flats([1 end],1:end)  = false;

if flag_nans
    % remove flat pixels bordering to nans
    flats(imdilate(log_nans,ones(3))) = false;
end

% identify sills
if nargout == 2;    
    % find sills and set marker
    Imr = -inf(size(dem));
    Imr(flats) = dem(flats);
    Imr = (imdilate(Imr,ones(3)) == dem) & ~flats;
    
    if flag_nans;
        Imr(log_nans) = false;
    end
end


