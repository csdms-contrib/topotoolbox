function cmap = ttscm(name,n)

%TTSCM scientific colormaps
%
% Syntax
%     
%     cmap = ttscm(name)
%     cmap = ttscm(name,n)
%     ttscm
%     allowedcmaps = ttscm
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
%     TTSCM without input arguments and one output argument returns the
%     names of the allowed colormaps as cell array.
%
% Input arguments
%
%     name           name of colormap (run ttscm without in- and output
%                    arguments for list)
%     n              number of colors in colormap    
%     
% Output arguments
%
%     cmap           n*3 colormap
%     allowedcmaps   allowed colormaps
%
% Example 1
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     C   = curvature(DEM);
%     lims = prcclip(C,2,true);
%     clr = ttscm('vik');
%     imageschs(DEM,C,'colormap',clr,'caxis',lims);
%
% Example 2: Show colormaps using the Big Tujunga Catchment data
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     cmaps = ttscm;
%     for r = 1:numel(cmaps); 
%        subplot(6,4,r); 
%        imageschs(DEM,[],'colormap', ttscm(cmaps{r}),...
%                         'colorbar',false,'ticklabels','none');
%        title(cmaps{r});
%     end
%
%     
% References: 
%
%      Crameri, F., (2018). Scientific colour-maps. Zenodo. 
%      http://doi.org/10.5281/zenodo.1243862
%
%      Crameri, F. (2018), Geodynamic diagnostics, scientific visualisation 
%      and StagLab 3.0, Geosci. Model Dev., 11, 2541-2562. 
%      doi:10.5194/gmd-11-2541-2018
%
%
% See also: IMAGESCHS, ttcmap
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 5. December, 2019

allowedcmaps = {'acton' 'bamako', 'batlow', 'hawaii', 'imola' 'nuuk' ...
                'buda' ...
                'devon','davos','oslo','bilbao','lajolla',...
                'grayC','broc','cork','vik','lisbon','tofino',...
				'berlin','turku','tokyo','lapaz','roma','oleron', ...
                'brocO','corkO','vikO'};

				
% get location of this function
p = fileparts(mfilename('fullpath'));
if nargin == 0 && nargout == 0
    figure('MenuBar','none','Color','w','Toolbar','none',...
           'Name', 'Scientific Colour Maps', 'NumberTitle', 'off');
           
	imshow([p filesep 'private' ...
            filesep '+ScientificColourMaps_FabioCrameri.png'],...
            'InitialMagnification','fit')
	return
elseif nargin == 0 && nargout == 1
    cmap = allowedcmaps;
    return
end

cmaptype = validatestring(name,allowedcmaps);
if nargin == 1
    n = 255;
else
    validateattributes(n,{'numeric'},{'>',1},'ttscm','n',2)
end

cmap = loadcmap(cmaptype);
ncolors = size(cmap,1);
if n ~= ncolors
	cmap    = interp1((1:ncolors)',cmap,linspace(1,ncolors,n));
end



