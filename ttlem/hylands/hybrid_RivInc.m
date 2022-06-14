
function [bedrock,sediment,lvSedTot,tot_DSY,Qs] =hybrid_RivInc(...
    bedrock, sediment,Lakes,DA,R,...
    dx2,K,K_bed_L,K_sed_L,m_area,n_slope,...
    i,k,dx_ik,boundaryNodesAll,mod_Domain,... 
    dt,phi,Ff,V,V_Lakes,H_star,pureIncision,...
    omega_cr,omega_cs)

% function part of HyLands: Hybrid Landscape evolution model, in
% this function the SPACE algorithm is applied and fluvial incision
% dynamics are being calculated. 
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



Q=DA.*R;
lvSedTot=0;
tot_DSY=0;
el=sediment+bedrock;
f_el=fillsinks(el);
f_el=f_el.Z;
el=el.Z;
lakes=f_el>el;
localD_i=lakes(i);
V_settling=ones(size(el)).*V;
V_settling(lakes(:))=V_Lakes;
a_r=K.*R.^m_area.*DA.^m_area;
a_r(isnan(el))=0;
a_r(a_r<0)=0;


if n_slope == 1
    Z=el;
    for r =numel(i):-1:1
        tt = a_r(i(r))*dt/dx_ik(r);
        Z(i(r)) = Z(i(r))-max((Z(i(r))-(Z(i(r)) + Z(k(r))*tt)./(1+tt)),0);
    end
else
    Z = funerosion_impnlin(el,dt_riv, a_r, uint64(i),uint64(k),dx_ik,n_slope);
end
ero=max(el-Z,0);
ero(Lakes.Z(:))=0;
if ~pureIncision
    
    omega_b=(K_bed_L./K).*ero;
    if omega_cr>0
        ero_bed=max(0,omega_b-omega_cr*dt*(1-exp(-omega_b./omega_cr*dt)));
    else
        ero_bed=omega_b;
    end
    omega_s=(K_sed_L./K).*ero;
    if omega_cs>0
        E_s_h=max(0,omega_s-omega_cs*dt*(1-exp(-omega_s./omega_cs*dt)));
    else
        E_s_h=omega_s;
    end
    
    
    E_bed=ero_bed.*(exp(-sediment.Z./H_star));
    clear ero_bed
    E_sed=E_s_h.*(1-exp(-sediment.Z./H_star));
    totEro=E_sed+E_bed;
    cells=(bedrock.Z(:)-totEro(:))<0;
    if sum(cells)
        scale=ones(size(totEro));
        scale(cells)=bedrock.Z(cells)./totEro(cells);
        clear totEro
        ero=ero.*scale;
        clear scale
        ero_bed=(K_bed_L./K).*ero;
        E_s_h=(K_sed_L./K).*ero;
        E_bed=ero_bed.*(exp(-sediment.Z./H_star));
        clear ero_bed
        E_sed=E_s_h.*(1-exp(-sediment.Z./H_star));
    end
    
    
    %% Calculate Qs
    % Use mex file to increase performance
%     [Qs,Qs_in]=UpdateSedFlux_mex(i,k,phi,E_sed,dx2,Ff,E_bed,V_settling,Q);
    [Qs,Qs_in]=UpdateSedFlux(i,k,phi,E_sed,dx2,Ff,E_bed,V_settling,Q);
    
    
    %% Update sediment thickness
    Ds=Qs./Q;
    Ds(i)=Ds(i).*V_settling(i);
    Ds(boundaryNodesAll)=0;
    
    loc500Plus=sediment.Z>500;
    iniThickness=sediment.Z(loc500Plus);
    sediment.Z(loc500Plus)=400;
    sed_T=H_star.*log((1./((Ds/(1-phi))./E_s_h-1)).*((exp((Ds/(1-phi)-E_s_h)/H_star)).*((((Ds/(1-phi))./E_s_h)-1).*exp(sediment.Z/H_star)+1)-1));
    
    sed_T(isnan(sed_T)|isinf(sed_T))=H_star*log(E_s_h(isnan(sed_T)|isinf(sed_T))./H_star + exp(sediment.Z(isnan(sed_T)|isinf(sed_T))/H_star));
    sed_T(E_s_h(:)<=0)=sediment.Z(E_s_h(:)<=0)+Ds(E_s_h(:)<=0)/(1-phi);
    sed_T(i(localD_i))=sediment.Z(i(localD_i))+Ds(i(localD_i))/(1-phi);
    clear Ds
    
    dHdt=sed_T-sediment.Z;
    sediment.Z(loc500Plus)=iniThickness;
    clear sed_T
    leavingSed=sum(Qs_in(:));
    sediment.Z=sediment.Z+dHdt;
    clear dHdt
    bedrock.Z=max(0,bedrock.Z-E_bed);
    lvSedTot=lvSedTot+leavingSed;
    tot_DSY=tot_DSY+sum(E_bed(mod_Domain)*Ff)*dx2;
    
    
else
    bedrock.Z=max(0,bedrock.Z-ero);
    lvSedTot=0;
    tot_DSY=0;
    Qs=0;
end

Qs=Qs/dt; %to m3/year

