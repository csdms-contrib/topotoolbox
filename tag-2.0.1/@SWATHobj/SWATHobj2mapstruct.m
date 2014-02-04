function [MS] = SWATHobj2mapstruct(SW,type)
% convert instance of SWATHobj to mappstruct
%
% Syntax
%
%     MS = SWATHobj2mapstruct(SW,type)
%
% Description
%
%     A mapstruct is a structure array that contains vector geographic
%     data. SWATHobj2mapstruct converts an instance of SWATHobj to a
%     mapstruct MS. MS can be exported to a shapefile using the function
%     shapewrite available with the Mapping Toolbox. It can be plotted with
%     the function mapshow.
%
%
% Input arguments
%
%     SW       SWATHobj
%     type     either 'lines' or 'points' can be exported from a SWATHobj
%
%     
% Output arguments
%
%     MS      mapstruct
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     SW = SWATHobj(DEM) % interactively create swath profile
%     [MS] = SWATHobj2mapstruct(SW,'lines');
%     figure, imageschs(DEM), hold on
%     mapshow(MS)
%
%
% See also: shapewrite
%
%
% Author: Dirk Scherler (scherler[at]caltech.edu)
% Date: 12. Dec, 2014



if ~isa(SW,'SWATHobj'); error('First input needs to be SWATHobj'); end


MS = struct;

switch type
    case 'points'
        
        ct = 1;
        for i = 1 : length(SW.Z)
            x = SW.X{i}; x = reshape(x,numel(x),1);
            y = SW.Y{i}; y = reshape(y,numel(y),1);
            for k = 1 : length(x)
                MS(ct).Geometry = 'Point';
                MS(ct).X = x(k);
                MS(ct).Y = y(k);
                ct = ct+1;
            end
        end
        
    case 'lines'
        
        ct = 1;
        for i = 1 : length(SW.Z)
            IM = true(size(SW.X{i}));
            B = bwboundaries(IM);
            ix = sub2ind(size(IM),B{1}(:,1),B{1}(:,2));
            x = SW.X{i}(ix);
            y = SW.Y{i}(ix);
            z = SW.Z{i}(ix);
            x = x(~isnan(z));
            y = y(~isnan(z));
            MS(ct).Geometry = 'Line';
            MS(ct).Type = 'Outline';
            MS(ct).X = x;
            MS(ct).Y = y;
            
            ct = ct+1;
            x = SW.xy0{i}(:,1);
            y = SW.xy0{i}(:,2);
            MS(ct).Geometry = 'Line';
            MS(ct).Type = 'Center';
            MS(ct).X = x;
            MS(ct).Y = y;
            
            ct = ct+1;
            x = SW.xy{i}(:,1);
            y = SW.xy{i}(:,2);
            MS(ct).Geometry = 'Line';
            MS(ct).Type = 'CenterFilt';
            MS(ct).X = x;
            MS(ct).Y = y;
            
            
        end
        
    otherwise
        error('Second input needs to be either ''lines'' or ''points''');
            
end











