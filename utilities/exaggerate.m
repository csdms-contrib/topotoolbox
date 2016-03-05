function exaggerate(axes_handle,exagfactor)

% elevation exaggeration in a 3D surface plot
%
% Syntax
%
%     exaggerate(axes_handle,exagfactor)
%
% Description
%
%     exaggerate is a simple wrapper for calling set(gca...). It controls
%     the data aspect ratio in a 3D plot and enables elevation
%     exaggeration.   
%
% Input
%
%     axes_handle   digital elevation model
%     exagfactor    exaggeration factor (default = 1)
%
% Example
%
%     load exampledem
%     for r = 1:4;
%     subplot(2,2,r);
%     surf(X,Y,dem); exaggerate(gca,r);
%     title(['exaggeration factor = ' num2str(r)]);
%     end
%
% 
% See also: SURF
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: June 2013

if nargin == 1;
    exagfactor = 1;
end
axis(axes_handle,'image');
set(axes_handle,'DataAspectRatio',[1 1 1/exagfactor]);


