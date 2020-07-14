function [bed, sed,Qs_Hill_OutD] = ...
    D_Hill(bed, sed,L_hill_i,mod_Domain,...
    i,k,fract,dx_ik,dt,dx,dx2,Qs_out_input,...
    phi,min_HillGrad)

% function part of HyLands: Hybrid Landscape evolution model 
% In this function deposition of hillslope material is caluclated following
% Carretier et al. 2016
%
% See also: HYLANDS_set, HYLANDS 
%
% =========================================================================
% Papers to cite when using HyLands or D_Hill function:
%
% * HyLands: Campforts B., Shobe M.C., et al. : HyLands 1.0: a Hybrid
% Landscape evolution model to simulate the impact of landslides and
% landslide-derived sediment on landscape evolution. Discussion paper in
% Geoscientific Model Development,
% https://geoscientific-model-development.net
%
% * Carretier, S., Martinod, P., Reich, M., & Godderis, Y. (2016).
% Modelling sediment clasts transport during landscape evolution. Earth
% Surface Dynamics, 4(1), 237–251. https://doi.org/10.5194/esurf-4-237-2016
%
% =========================================================================
%
% Author:   Benjamin Campforts (benjamin.campforts@gfz-potsdam.de)
%
% Date:     15. March, 2020


Qs_in=zeros(size(bed)); 
Qs_out=Qs_out_input*dt; % Flux in m3 per timestep
Qs_rem=Qs_out;
dH_Hill=zeros(size(bed));
H_i_temp=bed+ sed;
max_D=zeros(size(bed));

% No deposition allowed outside the modelled domain
max_dH=inf(size(bed));
max_dH(~mod_Domain)=0;

% Remaining Qs: for example when modelled domain is a non squarred polygon
%Qs_remaining
% totMassB=sum(sed(:))*dx2+sum(Qs_out(:))
for index = 1:numel(i) 
    %Of all sediment that is brought into the system Ff_Hill fraction of it
    %is immediatly brought into suspension and evacuated from the system. 
    dH=min(max_dH(i(index)),max(0,min(((Qs_in(i(index))/dx)/L_hill_i(index))/(1-phi),max_D(i(index)))));
    Qs_out(i(index))=Qs_out(i(index))+Qs_in(i(index))-dH*dx2*(1-phi);
    Qs_rem(i(index))=Qs_rem(i(index))+Qs_in(i(index))-dH*dx2*(1-phi);
    
    Qs_in(i(index))=0;% Next time we hit cell i, no more deposition
    
    Qs_in(k(index))=Qs_in(k(index))+fract(index)*Qs_out(i(index));    
    Qs_rem(i(index))=Qs_rem(i(index))-fract(index)*Qs_out(i(index));
    dH_Hill(i(index))  = dH_Hill(i(index))+dH;
    H_i_temp(i(index))=H_i_temp(i(index))+dH;    
    maxEl=max(H_i_temp(i(index))-min_HillGrad*double(dx_ik(index)),H_i_temp(k(index)));    
    max_D(k(index))=max(max_D(k(index)),max(0,maxEl-H_i_temp(k(index))));  
end

sed=sed+dH_Hill;

Qs_Hill_OutD=(sum(Qs_in(:))+sum(Qs_rem(:)))/dt;
% totMassA=sum(sed(:))*dx2+sum(Qs_Hill_OutD(:)*dt);diff=totMassB-totMassA
end
