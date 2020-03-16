function domainBorders=getDomainBorders(A,extend)

% function part of HyLands: Hybrid Landscape evolution model 
% If extend is true, the corders of the extend boundary domain are derived from the
% corners of the domain
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
% * SPACE: Shobe, C. M., Tucker, G. E., & Barnhart, K. R. (2017). The SPACE 1.0
% model: a Landlab component for 2-D calculation of sediment transport,
% bedrock erosion, and landscape evolution. Geoscientific Model
% Development, 10(12), 4577–4604. https://doi.org/10.5194/gmd-10-4577-2017
%
% * TTLEM: Campforts, B., Schwanghart, W., & Govers, G. (2017). 
% Accurate simulation of transient landscape evolution
% by eliminating numerical diffusion: the TTLEM 1.0 model. Earth Surface
% Dynamics, 5(1), 47–66. https://doi.org/10.5194/esurf-5-47-2017
%
% =========================================================================
%
% Author:   Benjamin Campforts (benjamin.campforts@gfz-potsdam.de)
%
% Date:     15. March, 2020


if extend
    domainBorders=[...
        [(size(A,1)+2) (size(A,1)+(2:size(A,1)-1)) (2*size(A,1)-1) ] ...                                     %LEFT
        [(size(A,1)*(size(A,2)-2)+2) (size(A,1)*(size(A,2)-2)+2:size(A,1)*(size(A,2)-1)-1) (size(A,1)*(size(A,2)-1)-1)] ...   %RIGHT
        2+size(A,1):size(A,1):size(A,1)*(size(A,2)-2)+2 ... %TOP
        (1+1*size(A,1):size(A,1):size(A,1)*(size(A,2)-2)+2)+size(A,1)-2];             %BOTTOM
    
else
    domainBorders=[...
        size(A,1)+(2:size(A,1)-1) ...                                     %LEFT
        (size(A,1)*(size(A,2)-1)+2:size(A,1)*size(A,2)-1)-size(A,1) ...   %RIGHT
        2+2*size(A,1):size(A,1):size(A,1)*(size(A,2)-3)+2 ... %TOP
        -1+size(A,1)*3:size(A,1):size(A,1)*(size(A,2)-2)-1];             %BOTTOM
end
end