function MASK = createmask(DEM,usehillshade)

%CREATEMASK create a binary mask using polygon mapping
%
% Syntax
%
%     MASK = createmask(DEM)
%     MASK = createmask(DEM,usehillshade)
%
% Description
%
%     createmask is an interactive tool to create a mask based on an
%     interactively mapped polygon.
%
% Input arguments
%
%     DEM           GRIDobj
%     usehillshade  use hillshade as background image ({false} or true)
%
% Output arguments
%
%     MASK   GRIDobj with logical mask
%
%
% See also: imroi, GRIDobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 17. August, 2017

if nargin == 1
    usehillshade = false;
end

figure
if usehillshade
   imageschs(DEM);
else
   imagesc(DEM);
end
title('Create polygon');
usedrawpolygon = ~verLessThan('matlab','9.5');
if usedrawpolygon
    ext = getextent(DEM);  

    h = drawpolygon('DrawingArea',...
        [ext(1) ext(3) ext(2)-ext(1) ext(4)-ext(3)]);
    pos = customWait(h);
else
    h = impoly; 
    pos = wait(h);
end

MASK = DEM;
MASK.name = 'mask';
MASK.Z = createMask(h);

if usedrawpolygon
    pos = h.Position;
else
    pos = getPosition(h);
end
delete(h);
hold on
plot(pos([1:end 1],1),pos([1:end 1],2));
hold off
title('done');


end


function pos = customWait(hROI)

% Listen for mouse clicks on the ROI
l = addlistener(hROI,'ROIClicked',@clickCallback);

% Block program execution
uiwait;

% Remove listener
delete(l);

% Return the current position
pos = hROI.Position;

end

function clickCallback(~,evt)

if strcmp(evt.SelectionType,'double')
    uiresume;
end

end
