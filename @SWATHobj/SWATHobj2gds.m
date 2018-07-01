function [MS] = SWATHobj2gds(SW,type)
%SWATHOBJ2GDS create a geographic data structure from a SWATHobj
%
% Syntax
%
%    MS = SWATHobj2gds(SW,type)
%
%
% Description
%
%     SWATHobj2gds creates a geographic data structure from a SWATHobj.
%     The geographic data structure can be exported as a shapefile. The
%     contents of the data structure are either points or lines, as
%     specified by the variable 'type'. If 'points' are selected, then the
%     position of all the elements of the matrix Z in the SWATHobj will be 
%     exported. If 'lines' are selected, then three lines are included: the
%     line that was used to create the swath (called 'Center' in the data
%     structure), the smoothened line (called 'CenterFilt' in the data
%     structure), and the outline of the swath (called 'Outline' in the data
%     structure).
%
%
% Input arguments
%
%     SW     Swath profile object (Class: SWATHobj)
%     type   string specifying the type of vector feature to include in the
%            geographic data structure ('points' or 'lines')
%
% Output
%
%     MS     geographic data structure
%
%
% Example 1
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     SW = SWATHobj(DEM);
%     MSp = SWATHobj2gds(SW,'points');
%     MSl = SWATHobj2gds(SW,'lines');
%     shapewrite(MSp,'swathobj_points.shp')
%     shapewrite(MSl,'swathobj_lines.shp')
%
%
% See also: SWATHobj
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: May, 2015



if ~isa(SW,'SWATHobj'); 
    error('First input needs to be SWATHobj'); 
end


MS = struct;

switch type
    case 'points'
        
        ct = 1;
        x = SW.X; x = reshape(x,numel(x),1);
        y = SW.Y; y = reshape(y,numel(y),1);
        for k = 1 : length(x)
            MS(ct).Geometry = 'Point';
            MS(ct).X = x(k);
            MS(ct).Y = y(k);
            ct = ct+1;
        end
        
    case 'lines'
        
        ct = 1;
        IM = true(size(SW.X));
        B = bwboundaries(IM);
        ix = sub2ind(size(IM),B{1}(:,1),B{1}(:,2));
        x = SW.X(ix);
        y = SW.Y(ix);
        z = SW.Z(ix);
        x = x(~isnan(z));
        y = y(~isnan(z));
        MS(ct).Geometry = 'Line';
        MS(ct).Type = 'Outline';
        MS(ct).X = x;
        MS(ct).Y = y;
        
        ct = ct+1;
        x = SW.xy0(:,1);
        y = SW.xy0(:,2);
        MS(ct).Geometry = 'Line';
        MS(ct).Type = 'Center';
        MS(ct).X = x;
        MS(ct).Y = y;
        
        ct = ct+1;
        x = SW.xy(:,1);
        y = SW.xy(:,2);
        MS(ct).Geometry = 'Line';
        MS(ct).Type = 'CenterFilt';
        MS(ct).X = x;
        MS(ct).Y = y;
        
        
    otherwise
        error('Second input needs to be either ''lines'' or ''points''');
        
end











