function boundaryNodes=getBoundaryNodes(A)
% function part of HyLands: Hybrid Landscape evolution model 
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
%
% Author:   Benjamin Campforts (benjamin.campforts@gfz-potsdam.de)
%
% Date:     15. March, 2020

boundaryNodes=[...
    1:size(A,1) ... %LEFT
    size(A,1)*(size(A,2)-1)+1:size(A,1)*size(A,2) ... %RIGHT
    size(A,1)+1:size(A,1):size(A,1)*(size(A,2)-2)+1 ... %TOP
    size(A,1)*2:size(A,1):size(A,1)*(size(A,2)-1)]; %BOTTOM
end