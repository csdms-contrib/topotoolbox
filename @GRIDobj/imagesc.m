function h = imagesc(DEM,varargin)

%IMAGESC Scale data in GRIDobj and display as image object
%
% Syntax
%
%     h = imagesc(DEM,varargin)
%
% Description
%
%     see imagesc 
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 18. February, 2022

%% Update on 18. Feb 2022
% nans in the data are set to transparent

[x,y] = refmat2XY(DEM.refmat,DEM.size);
ht = imagesc(x,y,DEM.Z,varargin{:});


axis xy
axis image

ht.AlphaData = ~isnan(DEM.Z);

if nargout == 1;
    h = ht;
end
