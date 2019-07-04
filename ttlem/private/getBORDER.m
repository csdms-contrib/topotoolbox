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
    case 'bl'
        BORDER.Z(end,:) = true;
        BORDER.Z(:,1) = true;
    case 'br'
        BORDER.Z(end,:) = true;
        BORDER.Z(:,end) = true;
    case 'tb'
        BORDER.Z(1,:) = true;
        BORDER.Z(end,:) = true;
    case 'tl'
        BORDER.Z(1,:) = true;
        BORDER.Z(:,1) = true;
    case 'tr'
        BORDER.Z(1,:) = true;
        BORDER.Z(:,end) = true;
    case 'lr'
        BORDER.Z(:,1) = true;
        BORDER.Z(:,end) = true;
    case 'tbl'
        BORDER.Z(1,:) = true;
        BORDER.Z(end,:) = true;
        BORDER.Z(:,1) = true;
    case 'tbr'
        BORDER.Z(1,:) = true;
        BORDER.Z(end,:) = true;
        BORDER.Z(:,end) = true;
        
    case 'tblr'
        BORDER.Z(:,1) = true;
        BORDER.Z(:,end) = true;
        BORDER.Z(1,:) = true;
        BORDER.Z(end,:) = true;
end
BORDER.Z = BORDER.Z*10000;
end

