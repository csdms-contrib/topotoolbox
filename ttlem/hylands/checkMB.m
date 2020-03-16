function MB=checkMB(iter,bed_Before,sed_Before,upl,dt,bedrock, sediment, ...
    MB_R,MB_H,mod_Domain,phi,MB_lim,Err_Message,dispMB)

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



if isempty(MB_R)
    MB_R=0;
end
if isempty(MB_H)
    MB_H=0;
end
if ~isscalar(upl)
    upl=upl(mod_Domain);
end


Diff_bed=sum(bed_Before-bedrock(mod_Domain));
Diff_reg=sum(sed_Before-sediment(mod_Domain))*(1-phi);

MB=Diff_bed+Diff_reg+sum(dt.*upl)-MB_R-MB_H;
if abs(MB)>MB_lim
    disp(['Iter equals ' num2str(iter)])
    disp(['MB equals ' num2str(MB)])
    error(['MB error after: ' Err_Message]);
end




if dispMB
    disp(['MB after ' Err_Message ' equals: ' num2str(MB)]);
end
if any(bedrock(:)<0)
    disp(['Iter: ' num2str(iter)]);
    error(['Negative values for bedrock after:' Err_Message]);
end
if any(sediment(:)<0)
    disp(['Iter: ' num2str(iter)]);
    error(['Negative values for sediment after:' Err_Message]);
end
if any(isnan(bedrock(:)))||any(isinf(bedrock(:)))
    disp(['Iter: ' num2str(iter)]);
    error(['Nan or Inf values for bedrock after:' Err_Message]);
end
if any(isnan(sediment(:)))||any(isinf(sediment(:)))
    disp(['Iter: ' num2str(iter)]);
    error(['Nan or Inf values for sediment after:' Err_Message]);
end


