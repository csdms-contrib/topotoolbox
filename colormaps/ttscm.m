function cmap = ttscm(name,varargin)

%TTSCM scientific colormaps
%
% Syntax
%     
%     cmap = ttscm(name)
%     cmap = ttscm(name,n)
%     ttscm
%
% Description
%
%     TTSCM provides access to scientific colormaps compiled by Fabio 
%     Crameri (http://www.fabiocrameri.ch/colourmaps.php). 
%
%     TTSCM(name) returns a 255*3 matrix with rgb values.
%
%     TTSCM(name,n) returns a n*3 matrix with rgb values.
%     
%     TTSCM without input and output arguments shows a figure with
%     available colormaps.
%
% Input arguments
%
%     zlimits        two element vector with maximum and minimum elevation
%     DEM            GRIDobj from which zlimits will be calculated    
%     
% Output arguments
%
%     cmap      n*3 colormap
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     C   = curvature(DEM);
%     lims = prcclip(C,2,true);
%     clr = ttscm('vik');
%     imageschs(DEM,C,'colormap',clr,'caxis',lims);
%     
% References: 
%
%      Crameri, F., (2018). Scientific colour-maps. Zenodo. 
%      http://doi.org/10.5281/zenodo.1243862
%
%      Crameri, F., Geodynamic diagnostics, scientific visualisation 
%      and StagLab 3.0, Geosci. Model Dev. Discuss., doi:10.5194/gmd-2017-328, 
%      in open review, 2018.
%
%
% See also: IMAGESCHS, ttcmap
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 14. June, 2018

allowedcmaps = {'devon','davos','oslo','bilbao','lajolla',...
                'grayC','broc','cork','vik','lisbon','tofino',...
				'berlin','turku','tokyo','lapaz','roma','oleron'};

				
% get location of this function
p = fileparts(mfilename('fullpath'));
if nargin == 0
    figure('MenuBar','none','Color','w','Toolbar','none',...
           'Name', 'Scientific Colour Maps', 'NumberTitle', 'off');
           
	imshow([p filesep 'private' ...
            filesep '+ScientificColourMaps_FabioCrameri.png'],...
            'InitialMagnification','fit')
	if nargout == 1
	    cmap = [];
	end
	return
end

cmaptype = validatestring(name,allowedcmaps);
if nargin == 1
    n = 255;
else
    validateattributes(n,{'numeric'},{'>',1},ttscm,'n',2)
end

cmap = loadcmap(cmaptype);
ncolors = size(cmap,1);
if n ~= ncolors
	cmap    = interp1((1:ncolors)',cmap,linspace(1,ncolors,n));
end



