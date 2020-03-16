function [sediment, bedrock,Qs_change_LS,suspended_Sed,...
    LS_properties,slidePlain,LSloc,LS_Total_E] = LS_cullman(...
    sediment, bedrock,nb_M,...
    locationsD,X,Y,phi_Sc,mod_Domain_M,slidePlain,FD_single,dx_ik,...
    dt,dx2,iter,maxLS_Size,phi,Ff_LS,C_eff,g,rho,LS_DtChar,...
    save_LS_Data,resultsdir, fileprefix,LS_Total_E,verbose,saveeach) 

% Function part of HyLands: Hybrid Landscape evolution model: in this
% function, the location of critical landslide nodes is identified based on
% there probability for sliding, as well as the magnitude of the landslide.
% The method is based on the Cullman appraoch, following Densmore et al.
% 1998. See referenced papers for further details. See HYLANDS and
% HYLANDS_set for info on paramter values.
%
% See also: HYLANDS, HYLANDS_set
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
% Other relevant references:
%
% * Carretier, S., Martinod, P., Reich, M., & Godderis, Y. (2016).
% Modelling sediment clasts transport during landscape evolution. Earth
% Surface Dynamics, 4(1), 237–251. https://doi.org/10.5194/esurf-4-237-2016
%
% * Densmore, A. L., Ellis, M. A., & Anderson, R. S. (1998). Landsliding
% and the evolution of normal-fault-bounded mountains. Journal of
% Geophysical Research: Solid Earth, 103(B7), 15203–15219.
% https://doi.org/10.1029/98JB00510
%
% =========================================================================
%
% Author:   Benjamin Campforts (benjamin.campforts@gfz-potsdam.de)
%
% Date:     15. March, 2020


%% Bedrock lS
% Define maximum landslide size
ElBefore=sediment+bedrock;
maxPixels=maxLS_Size/dx2;
el=bedrock.Z+sediment.Z;
phi_Sc=phi_Sc.Z;
H=el(FD_single.ix)-el(FD_single.ixc);
beta=H./dx_ik;
beta_d=atand(beta);
phi_d=atand(phi_Sc(FD_single.ix));
if ~isscalar(C_eff)
    C_eff=C_eff(FD_single.ix);
end
Hc=(4*C_eff/(g*rho)).*((sind(beta_d).*cosd(phi_d)) ./ (1-cosd(beta_d-phi_d)) );
Hc(beta_d<=phi_d)=Inf;
p=H./Hc;
p(beta_d<phi_d)=0; % No slides if slope is lower than threshold slope
slides=rand(size(p))<p;
i_slide= unique(FD_single.ixc(slides));%Critical nodes must not be unique, select unique ones only
i_slide=datasample(i_slide,round(numel(i_slide)*dt/LS_DtChar),'Replace',false);
Splain=zeros(size(bedrock.Z));
Splain(i_slide)=1;
Splain(~locationsD)=0;
i_slide=find(Splain);
if verbose
    disp(['Nb Deep-Seated bedrock LS = ' num2str(numel(i_slide))]);
end
clear H Hc
Qs_change_LS=zeros(size(el));
suspended_Sed=0;
LS_properties.Size=[];
LS_properties.Volume=[];
LS_properties.Volume_Bed=[];
LS_properties.Volume_Reg=[];

if save_LS_Data
    LS_Vol_Cur=zeros(size(slidePlain));
end

LSloc=i_slide;

while ~isempty(i_slide)
    ind=1;
    cP=i_slide(ind);
    cP_el=el(cP);
    nb=nb_M(cP,:);
    nb(nb==0)=[];
    nb_up=nb(el(nb)>cP_el);
    nb_up(mod_Domain_M(nb_up)==0)=[];
    distToIni_all= sqrt((X(cP)-X(nb_up)).^2 + (Y(cP)-Y(nb_up)).^2);
    all_iP_el=el(nb_up);
    s_slide_all=((max(all_iP_el-cP_el,0))./distToIni_all);
    phi_cp=phi_Sc(i_slide(ind));
    nb_up(s_slide_all<phi_cp)=[];
    s_slide_all(s_slide_all<phi_cp)=[];
    s_slide=double(max(s_slide_all));
    storeV_bed=0;
    storeV_sed=0;
    upstream=0;
    uP=nb_up';
    ang_sl=(phi_cp+s_slide)/2;
    
    while ~isempty(uP)&& upstream<=maxPixels
        distToIni= double(sqrt((X(cP)-X(uP(1))).^2 + (Y(cP)-Y(uP(1))).^2));
        newEl=max(cP_el,cP_el+distToIni*ang_sl);
        if newEl<el(uP(1))
            upstream=upstream+1;
            sed_LS_E=min(sediment.Z(uP(1)),el(uP(1))-newEl);
            sediment.Z(uP(1))=sediment.Z(uP(1))-sed_LS_E;
            bedrock.Z(uP(1))=newEl-sediment.Z(uP(1));
            storeV_sed=storeV_sed+sed_LS_E*(1-phi)*dx2;
            storeV_bed=storeV_bed+((el(uP(1))-newEl)-sed_LS_E)*dx2;
            el(uP(1))=newEl;
            nb_up=nb_M(uP(1),:);
            nb_up(nb_up==0)=[];
            nb_up(mod_Domain_M(nb_up)==0)=[];
            nb_up(el(nb_up)<el(uP(1)))=[];
            uP_l=nb_up;
            uP=[uP; uP_l'];
            [~,idx] = unique(uP,'first');
            uP = uP(sort(idx));
            indslide=(i_slide==uP(1));
            i_slide(indslide)=[];
            slidePlain(uP(1))=iter;
        end
        uP(1)=[];
    end
    storeV=storeV_sed+storeV_bed;
    Qs_change_LS(cP)=Qs_change_LS(cP)+(storeV)/dt*(1-Ff_LS);
    suspended_Sed=suspended_Sed+storeV*Ff_LS;
    LS_properties.Size=[LS_properties.Size upstream];
    LS_properties.Volume=[LS_properties.Volume storeV];
    LS_properties.Volume_Reg=[LS_properties.Volume_Reg storeV_sed];
    LS_properties.Volume_Bed=[LS_properties.Volume_Bed storeV_bed];
    if ~isempty(i_slide)
        i_slide(1)=[];
    end
end
LS_act_E=uint16(ElBefore.Z-(sediment.Z+bedrock.Z));
LS_Total_E=LS_Total_E+LS_act_E;

if  mod(iter,saveeach)==0&& save_LS_Data
    save([resultsdir fileprefix '_res_' num2str(round(bedrock.cellsize)) '_LS_loc_' num2str(round(iter*dt)) '.mat'],'LSloc');
    save([resultsdir fileprefix '_res_' num2str(round(bedrock.cellsize)) '_LS_act_E_' num2str(round(iter*dt)) '.mat'],'LS_act_E');
end
if any(el(:)<0)
    disp('negativeV')
    error('An error occured during landsliding. The elevation cannot be negative.')
end
end