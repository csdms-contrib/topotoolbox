function EGM96 = egm96heights(DEM)

%EGM96HEIGHTS read and resample EGM96 geoid heights
%
% Syntax
%
%     EGM96 = egm96heights(DEM)
%     egm96heights
%
% Description
%
%     EGM96HEIGHTS returns the egm 96 geoid heights as GRIDobj that
%     is spatially aligned with the input GRIDobj DEM.
%
%     This function requires the mapping toolbox that includes the 
%     function egm96geoid.
%
% Input arguments
%
%     DEM      GRIDobj
%     
% Output arguments
%
%     EGM96    GRIDobj with geoid heights
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     EGM96 = egm96heights(DEM);
%     imageschs(DEM,EGM96)
%
%  
% See also: GRIDobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 8. May, 2018

[Z, refvec] = egm96geoid(1,[90 -90],[-180 180]);
R = refvecToGeoRasterReference(refvec,size(Z));
lat = linspace(R.LatitudeLimits(1),R.LatitudeLimits(2),size(Z,1));
lon = linspace(R.LongitudeLimits(1),R.LongitudeLimits(2),size(Z,2));
EGM96 = GRIDobj(lon,lat,Z);

if nargin == 0 && nargout == 0
    figure
    imagesc(EGM96);
    load coast
    hold on
    plot(long,lat,'k');
    hold off
    
elseif nargin > 0

        switch DEM.georef.SpatialRef.CoordinateSystemType
            case 'geographic'
                EGM96 = resample(EGM96,DEM);
            otherwise
                EGM96 = reproject2utm(EGM96,DEM);
        end

end
