function [ic,icd,W0] = crossflats(dem,type)

% cross flats of a digital elevation model
%
% Syntax
%     
%     [IXf,IXn,W0] = crossflats(dem,type)
%
% Description
%
%     crossflats is a subroutine used by some of the flowdirection 
%     algorithms in the toolbox. While routeflats recursively creates flow 
%     paths through flat terrain, crossflats calculates the connectivity
%     between nodes draining into flat nodes and nodes draining the flats. 
%     The advantage compared to routeflats is low computational costs 
%     in particular for very large flats (such as huge water bodies in a
%     DEM). The disadvantage of this method is, that, when creating a flow
%     direction matrix with the output of crossflats, the nodes in flat
%     terrain will be lacking in analysis such as flowpathbuffer or
%     flowpathanalysis.
%
%     The user may choose between single and multiple flow routing.
%     Multiple flow routing is not recommended in cases when the flats have
%     very many outlets.
%
% Input
%
%     dem       Digital elevation model
%     type      'single' (default) or 'multi' flowdirection
%
% Output
%     
%     IXf,IXn   vectors with linear indices of cells surrounding flats. IXf
%               are cells upstream and IXn are cells downstream to flat
%               terrain.
%     W0        W0 is a matrix same size as dem, where each element is one
%               or greater, when the cell is an outlet of a flat. Elements
%               greater than one refer to the number of cells in flats that
%               should be added to the outlet.
%
%
% See also: ROUTEFLATS
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 9. July, 2009


% check image processing toolbox version
if verLessThan('images','6.3')
    warning('TopoToolbox:version',...
          'crossflats works best with ITP versions equal or higher 6.3')
    flagver = 0;
else
    flagver = 1;
end


% check input arguments
if nargin == 1;
    type = 'single';
end
    

mindem = min(dem(:));
siz    = size(dem);
% border cells

nhood = ones(3,3);

% handle NaNs
log_nans = isnan(dem);
if any(log_nans(:));
    flag_nans = true;
else
    flag_nans = false;
end

% if there are nans.
if flag_nans
    fillval       = mindem-1;
    dem(log_nans) = fillval;
end

% identify flats
% flats: logical matrix with true where cells don't have lower neighbors
if flag_nans
    flats = imerode(dem,nhood) == dem & ~log_nans;
    flats(imdilate(log_nans,nhood))= false;
else
    flats = imerode(dem,nhood) == dem;
end

% remove flats located on the edge of the dem
flats(1:end,[1 end]) = false;
flats([1 end],1:end) = false;

if flagver
    % if you have the image processing toolbox version 6.3 or later,
    % crossflats uses bwconncomp. bwconncomp is much faster than the
    % workaround with bwlabel and regionprops. Still, in order to support
    % older versions, too, the bwlabel alternative is done for earlier IPT
    % versions.
    CC    = bwconncomp(flats);
else
    SS    = regionprops(bwlabel(flats),'PixelIdxList'); %#ok
    CC.PixelIdxList  = {SS.PixelIdxList};  
    clear SS
end

% perimeter pixel of flats
P = bwperim(flats,8);

% find cells in the perimeter pixel of flats that have 
% neighbors that are not flats and have equal or lower elevation.
shiftval = [1 -1 siz(1) -siz(1) ...
            1+siz(1) -1+siz(1) -1-siz(1) 1-siz(1)];
        
% this section deals with identifying the perimeter pixels of the flats and
% their linkage with the outlet cell of flats. 
switch type
    case 'single'
        % pixel values of perimeter in each flat
        [CC.PixelIdxList, CC.Area] = cellfun(@(IX) perimpixels(IX),CC.PixelIdxList',...
                                     'UniformOutput',false);
        
        [CC.OutletIdxList,CC.PixelIdxList] = ... 
            cellfun(@(IX) findoutletsingle(IX),CC.PixelIdxList,'UniformOutput',false);
        CC.Area = cellfun(@(area,IX) adjustarea(area,IX),CC.Area,CC.PixelIdxList,'UniformOutput',false);
    case 'multi'
        
        [CC.PixelIdxList,CC.Area] = cellfun(@(IX) perimpixels(IX),CC.PixelIdxList',...
                                     'UniformOutput',false);
        [CC.OutletIdxList,CC.PixelIdxList] = ... 
            cellfun(@(IX) findoutletmulti(IX),CC.PixelIdxList,'UniformOutput',false);
        CC.Area = cellfun(@(area,IX) adjustarea(area,IX),CC.Area,CC.PixelIdxList,'UniformOutput',false);
end

% create ic and icd vectors
ic  = cell2mat(CC.PixelIdxList);
icd = cell2mat(CC.OutletIdxList);
W0  = cell2mat(CC.Area);

% and this is it





% some subfunctions
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [IXoutlet,IXperim] = findoutletsingle(IXperim)
% findoutletsingle searches for all the perimeter pixels in a flat the
% outlet pixel. It than connects the outlet pixel with the perimeter pixels
%


r = 1;
IXoutlet = [];
while isempty(IXoutlet) && r <= 8;
    
    IXtest   = IXperim+shiftval(r);
    IXoutlet = IXtest(find((~(flats(IXtest)) & dem(IXtest)<=dem(IXperim)),1,'first')); 
    r = r+1;
end

if isempty(IXoutlet);
    IXperim = [];
else
    IXoutlet = repmat(IXoutlet,numel(IXperim),1);
end
end


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [IXoutlet,IXperim] = findoutletmulti(IXperim)
% same as findoutletsingle, except that it works for the multiple flow
% direction algorithm. The major difference between both is, that
% findoutletmulti determines a pixel where all perimeter pixels drain into
% and which then drains towards all outlet pixels. In case of only one
% outlet pixel all perimeter pixels only drain towards this outlet pixel.

IXoutlet = [];
for r = 1:8;
    
    IXtest   = IXperim+shiftval(r);
    IXoutlet = [IXoutlet;IXtest(~(flats(IXtest)) & dem(IXtest)<=dem(IXperim))]; %#ok

end 

if isempty(IXoutlet);
    IXperim = [];
else
    if numel(IXoutlet)== 1;
        IXoutlet = repmat(IXoutlet,numel(IXperim),1);
    else
        
        IXoutlet = [repmat(IXoutlet(1),numel(IXperim),1);...
                    IXoutlet(2:end)];
        IXperim  = [IXperim;repmat(IXoutlet(1),numel(IXoutlet)-numel(IXperim),1)];
    end
end
end
    


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [IXperim,area] = perimpixels(IX)
% perimpixels returns the pixel indices on the perimeter of a flat region
    
IXperim = IX(P(IX));
area    = numel(IX)-numel(IXperim);

end

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function W0 = adjustarea(area,IX)
% adjustarea calculates the drainage area for each perimeter pixel thereby
% incorporating the total area of the flat region.

if isempty(IX)
    W0 = [];
else
    nrIX = numel(IX);
    W0 = 1 + repmat(area./nrIX,nrIX,1);
end
end

end




