function [Qs,Qs_in]=UpdateSedFlux(i,k,phi,E_reg,dx2,Ff,E_bed,V,Q)
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

Qs=zeros(size(E_bed));
Qs_in=zeros(size(E_bed));
for ii = 1:numel(i)    
    Qs(i(ii))=(Qs_in(i(ii))+(1-phi)*E_reg(i(ii))*dx2+(1-Ff)*E_bed(i(ii))*dx2)/(1+V(i(ii))*dx2/(Q(i(ii))));
    Qs_in(k(ii))=Qs_in(k(ii))+Qs(i(ii));
    Qs_in(i(ii))=0;
end
end
