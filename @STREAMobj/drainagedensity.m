function d = drainagedensity(S,FD,outputformat)

%DRAINAGEDENSITY drainage density of a stream network
%
% Syntax
%
%     D = drainagedensity(S,FD)
%     d = drainagedensity(S,FD,outputformat)
%     d = drainagedensity(S,FD,'nal')
%     d = drainagedensity(S,FD,ix)
%
% Description
%
%     DRAINAGEDENSITY calculates the drainage density of a stream network
%     which is calculated by dividing the stream length by the upslope area
%     at each location.
%
%     D = drainagedensity(S,FD) returns a GRIDobj where non-stream pixels
%     are nan and stream pixels contain their upstream drainage density.
%     This is the same as drainagedensity(S,FD,'output','grid')
%
%     d = drainagedensity(S,FD,'nal') returns the drainage densities for
%     each stream location as a node-attribute list.
%
%     d = drainagedensity(S,FD,ix) takes the linear index ix into the DEM
%     from which S and FD were derived and returns the drainage densities
%     only for these locations.
%
% Input arguments
%
%     S             stream network (STREAMobj)
%     FD            flow directions (FLOWobj)
%     outputformat  'grid', 'nal' or linear index ix
%
% Output arguments
%
%     D or d        GRIDobj, nal or vector with drainage densities
%
% Example 1: Calculate drainage densities for selected catchments
%     
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     [DB,IX] = drainagebasins(FD,S);
%     d = drainagedensity(S,FD,IX);
%     DD = GRIDobj(DEM);
%     DD.Z(DB.Z~=0) = d(DB.Z(DB.Z~=0));
%     imageschs(DEM,DD)
%
% Example 2: Calculate grid with drainage densities
%
%     D = drainagedensity(S,FD)
%     imageschs(S,D)
%
% Example 3: Plot stream profile together with drainage densities
%
%     d = drainagedensity(S,FD,'nal');
%     plotdz(S,DEM,'color',d,'LineWidth',2);
%
% See also: STREAMobj, FLOWobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 6. September, 2017

if nargin == 2
    outputformat = 'grid';
end

a = flowacc(FD)*FD.cellsize^2;
a = getnal(S,a);
d = distance(S,'accumdownstream')./a;

if ischar(outputformat)
    outputformat = validatestring(outputformat,{'grid','nal'});
else
    outputlocations = outputformat;
    outputformat = 'ix'; 
end

switch outputformat
    case 'grid'
        d = STREAMobj2GRIDobj(S,d);
    case 'nal'
       
    case 'ix'
        [I,locb] = ismember(outputlocations,S.IXgrid);
        dd = nan(numel(outputlocations),1);
        dd(I) = d(locb(I));
        d  = dd;
        if any(~I)
            warning('Some locations are not on the channel network.')
        end
end
