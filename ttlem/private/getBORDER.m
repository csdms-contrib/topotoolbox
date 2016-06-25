function BORDER = getBORDER(DEM,p)
%Function to get the border in fuction of the boundary condtions. 
%
% Syntax
%
%       BORDER = getBORDER(DEM,p)
%
% Description
%
% Function to get the border in fuction of the boundary condtions.
%
% Input
%
%       DEM        DEM (digital elevation model) (GRIDobj)
%       p          structure array with parameter definitions (see ttlemset)
% Output
%
%       BORDER
%
% Example
%
%
% See also:
%
% Authors:  Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
%           Benjamin Campforts (benjamin.campforts@ees.kuleuven.be)
%
%
% Date: 28. Januari, 2015

BORDER = GRIDobj(DEM,'logical');
switch p.FlowBC
    case 'b'        
        BORDER.Z(end,:) = true;
    case 't'
        BORDER.Z(1,:) = true;
    case 'r'
        BORDER.Z(:,end) = true;
    case 'l'
        BORDER.Z(:,1) = true;
    case 'lb'
        BORDER.Z(end,:) = true;
        BORDER.Z(:,1) = true;
    case 'rb'
        BORDER.Z(end,:) = true;
        BORDER.Z(:,end) = true;
    case 'all'
        BORDER.Z(:,1) = true;
        BORDER.Z(:,end) = true;
        BORDER.Z(1,:) = true;
        BORDER.Z(end,:) = true;
end
BORDER.Z = BORDER.Z*10000;
end

