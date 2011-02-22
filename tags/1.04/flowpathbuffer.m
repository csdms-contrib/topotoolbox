function IX = flowpathbuffer(siz,ord,pixdist)

% create a buffer around a stream and extract indices of cells
%
% Syntax
%
%     IXc = flowpathbuffer(siz,IX,pixdist)
%
% Description
%
%     flowpathbuffer extracts linear indices of cells that lie within a
%     specified distance perpendicular to the flow path. This function is 
%     particularly useful, when analyzing the slopes adjacent to a channel.
%     The output of flowpathbuffer is a cell array same size as IX where
%     each cell contains the linear indices of the cells within the buffer.
%     See the example below and the function cellfun when you are new to
%     working with cell arrays.
%
% Input
%
%     siz       size of the digital elevation model as return by size(dem)
%     ord       linear index of stream as returned by flowpathextract
%     pixdist   buffer distance (unit: # of gridcells) 
%
% Output
%
%     IXc       cell array same size as ord containing the linear
%               indices of cells in the DEM (or whatever) that are
%               attributed to IX.
%
% Example
%
%     load exampleDEM
%     M  = flowdir_single(dem);
%     % find cell with highest elevation
%     [ignore,channelstart] = max(dem(:));
%     % calculate flowpath from highest point
%     IX = flowpathextract(M,channelstart);
%     % plot dem and flowpath location
%     imagesc(X(1,:),Y(:,2),dem); axis image; axis xy
%     hold on
%     plot(X(IX),Y(IX),'k')
%     hold off
% 
%     % calculate the distance downstream
%     dis    = [0; cumsum(hypot(X(IX(1:end-1))-X(IX(2:end)),...
%                  Y(IX(1:end-1))-Y(IX(2:end))))];
% 
%     % calculate mean slope in a distance of 3 pixels
%     % perpendicular to the flowpath
%     G       = gradient8(dem,abs(Y(1)-Y(2)));
%     IXslope = flowpathbuffer(size(dem),IX,3);
%     S       = cell2mat(cellfun(@(x) mean(G(x)),...
%               IXslope,'uniformoutput',false));
%
%     % plot results
%     figure
%     [AX,H1,H2] = plotyy(dis,dem(IX),dis,S);
%
%     % annotate the plot 
%     set(get(AX(1),'Ylabel'),'String','elevation') 
%     set(get(AX(2),'Ylabel'),'String','mean slope adjacent to channel')  
%     xlabel('distance downstream [m]')
%
%
% See also: FLOWPATHEXTRACT, CELLFUN, FLOWDISTANCE
% 
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 15. March, 2009


nrc     = prod(siz);
BW      = false(siz);
BW(ord) = true;
[D,L]   = bwdist(BW);

I       = D(:)>pixdist;
L       = L(:);
L(I)    = [];
IXc      = (1:nrc)';
IXc(I)   = [];

IX       = accumarray(L,IXc,[nrc,1],@(x) {x});

IX = IX(ord);

 