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
% See also: HYLANDS_set, HYLANDS
%
% =========================================================================
% Papers to cite when using HyLands:
%
% * HyLands: Campforts B., Shobe M.C., et al. : HyLands 1.0: a Hybrid
% Landscape evolution model to simulate the impact of landslides and
% landslide-derived sediment on landscape evolution. Discussion paper in
% Geoscientific Model Development,
% https://geoscientific-model-development.net
%
% * TTLEM: Campforts, B., Schwanghart, W., & Govers, G. (2017). 
% Accurate simulation of transient landscape evolution
% by eliminating numerical diffusion: the TTLEM 1.0 model. Earth Surface
% Dynamics, 5(1), 47–66. https://doi.org/10.5194/esurf-5-47-2017
%
% =========================================================================
% Authors:  Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
%           Benjamin Campforts (benjamin.campforts@gfz-potsdam.de)
%
% Date:     15. March, 2020

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
        
    case  'ul_cor'
        BORDER.Z(1,:) = true;
        BORDER.Z(end,:) = true;
        BORDER.Z(:,end) = true;
        BORDER.Z(:,1) = true;
        BORDER.Z(1,1) = false;
    case  'll_cor'
        BORDER.Z(1,:) = true;
        BORDER.Z(end,:) = true;
        BORDER.Z(:,end) = true;
        BORDER.Z(:,1) = true;
        BORDER.Z(end,1) = false;
    case  'ur_cor'
        BORDER.Z(1,:) = true;
        BORDER.Z(end,:) = true;
        BORDER.Z(:,end) = true;
        BORDER.Z(:,1) = true;
        BORDER.Z(1,end) = false;
    case  'lr_cor'
        BORDER.Z(1,:) = true;
        BORDER.Z(end,:) = true;
        BORDER.Z(:,end) = true;
        BORDER.Z(:,1) = true;
        BORDER.Z(end,end) = false;
end
corners={'ul_cor','ll_cor','ur_cor','lr_cor'};
if ~any(strcmp(corners,p.FlowBC))
%Check corners
if BORDER.Z(1,2)==0||BORDER.Z(2,1)==0
    BORDER.Z(1,1)=0;
end
if BORDER.Z(1,end-1)==0||BORDER.Z(2,end)==0
    BORDER.Z(1,end)=0;
end
if BORDER.Z(end-1,1)==0||BORDER.Z(end,2)==0
    BORDER.Z(end,1)=0;
end
if BORDER.Z(end-1,end)==0||BORDER.Z(end,end-1)==0
    BORDER.Z(end,end)=0;
end
end
BORDER.Z = BORDER.Z*10000;
end

