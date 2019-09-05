function geotiffextent(filename)

%GEOTIFFEXTENT plot extent of geotiffs in a folder
%
% Syntax
%
%     geotiffextent(folder)
%
% Description
%
%     This tiny tool takes a folder as input and returns a plot with
%     patches that indicate the extent of each geotiff in the folder. Note 
%     that currently the function only reads the bounding box of each tif.
%     Thus, in order for the function to properly work, all tifs in the
%     directory should have the same coordinate system.
%
% Input arguments
%
%     folder     folder with tifs. 
%
% Output arguments
%
%     none
%
%
% See also: GRIDobj, geotiffinfo
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. February, 2019

if isfolder(filename)
    files = dir([filename filesep '*.tif']);
    x  = cell(size(files));
    y  = cell(size(files));
    parfevalOnAll(gcp(), @warning, 0, 'off');
    parfor r = 1:numel(files)
       info = geotiffinfo([files(r).folder filesep files(r).name]);
       xx   = info.BoundingBox(:,1);
       yy   = info.BoundingBox(:,2);
       x{r} = [xx(1) xx(2) xx(2) xx(1) xx(1)]';
       y{r} = [yy(1) yy(1) yy(2) yy(2) yy(1)]';       
    end
    parfevalOnAll(gcp(), @warning, 0, 'on');
end

X = horzcat(x{:});
Y = horzcat(y{:});

patch(X,Y,[.9 .9 .9]);
    
