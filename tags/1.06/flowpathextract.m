function ord = flowpathextract(M,channelstart)

% extract linear indices of a single flowpath in a DEM
%
% Syntax
%
%     IX = flowpathextract(M,channelstart)
%
% Description
%
%     flowpathextract finds the linear indices of a flow path starting at
%     channelstart based on the single flow direction matrix M.
%
% Input
%
%     flowdir       sparse SINGLE flow direction matrix
%     channelstart  linear index (scalar) where channel starts
%
% Output
%
%     IX            linear index
%
% Example
% 
%     load exampleDEM
%     M = flowdir_single(dem);
%     % find cell with highest elevation
%     [ignore,channelstart] = max(dem(:));
%     IX = flowpathextract(M,channelstart);
%     imagesc(X(1,:),Y(:,2),dem); axis image; axis xy
%     hold on
%     plot(X(IX),Y(IX),'k')
%
%     % calculate the distance downstream
%     dis    = [0; cumsum(hypot(X(IX(1:end-1))-X(IX(2:end)),...
%                    Y(IX(1:end-1))-Y(IX(2:end))))];
%     figure
%     plot(dis,dem(IX))
%     xlabel('distance downstream [m]')
%     ylabel('elevation')
% 
%
% See also: FLOWDIR, GETPTS
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 15. March, 2009



% check if single flow direction Matrix
% is used
if any(sum(spones(M),2)>1);
    error('TopoToolbox:incorrectinput',...
          'single flow direction matrix must be used')
end


% nr of cells
nrc = size(M,1);


z = zeros(nrc,1);
z(channelstart)=1;
M = (speye(nrc)-M');
f = M\z;
f = M\f;
ord  = (1:nrc)';
[f,IX] = sort(f,'ascend');
ord  = ord(IX(f~=0));


