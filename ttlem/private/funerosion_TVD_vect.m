function Z= funerosion_TVD_vect(p,Z,dt,A,i,k,dx_ik,kk,ii,dx_centered,upl)
% Explicit solution for river incision, calucalted by solving the Stream Power
% Law.
%
% Syntax
%
%       Z = funerosion_ex(p,Z,dt,dx, A, i,k,dx_ik)
%
% Description
%
%       Explicit solution for river incision, calucalted by solving the
%       Stream Power Law. Internally, time steps are adapted as to honor
%       the cfl criterion.
%
% Input
%
%       p         parameter values
%       Z         digital elevation model (matrix)
%       dt        time step
%       dx        spatial step
%       A         velocity field
%       i         giver, derived from the flow
%                 direction FD.ix where ix is an edge attribute and
%                 represetns topologically sorted nodes (givers)
%       k         receiver, derivedfrom the flow direction FD.ixc where
%                 ixc is an edge attribute and represetns topologically
%                 sorted nodes (receivers) dx_ik horizontal distance
%       dx_ik     distance between giver and receiver 
%       kk        Receiver of receiver 
%       ii        Giver of giver 
%       dx_centered distance between giver of giver and receiver
%       upl       Uplift rate (m/yr)
%
% Output
%
%       Z       digital elevation model with adapted river elevation
%
% Example
%
%
% See also:
%
% Authors: Benjamin Campforts (benjamin.campforts[at]ees.kuleuven.be)
%          Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
%
% Date: 28. Januari, 2015



a_ori=-A(i);
a_m = min(0,a_ori);
a_p = max(0,a_ori);
a_vect=a_ori;

nrtsteps = ceil(dt/(p.cfls*min(dx_ik)/max(abs(a_vect))));
dte = dt/nrtsteps;
time=dt;

while time>0;
    %% TVD
    % Nonlinear solution, calculated by integrating the explicit FDM
    % solution to the incision speed
    if p.n~=1
        s1=((max(Z(i)-Z(k),0))./dx_ik); %dx_ik? 
        %Very low slopes are leading through very small timesteps so they are set to 0.
        s1(s1<1e-4)=0;
        exp_f=s1.^(p.n-1);
        exp_f(isinf(exp_f))=1;
        a_vect(exp_f~=0)=a_ori(exp_f~=0).*exp_f(exp_f~=0);
        a_m = min(0,a_vect);
        a_p = max(0,a_vect);
        
        % Check timestep
        dt_calc = p.cfls*min(dx_ik)/max(abs(a_vect));
        if dt_calc<dte
            disp(['For stability, dte is set to: ' num2str(dt_calc)]);
            dte=dt_calc;
        end
    end
    %Time
    time=time-dte;
    if time<0
        dte=dte+time;
        time=0;
    end
       
   
        el_c=Z(i);
        el_d=Z(k);
        kk_nb=kk;
            kk_nb(isnan(kk))=1;
            el_d2=Z(kk_nb);    
            el_d2(isnan(kk))=nan;%el_d(isnan(kk));
            iiNb=ii;
            iiNb(isnan(ii))=1;
            el_up=Z(iiNb);
            el_up(isnan(ii))=nan;
        r_TVD=(el_d2-el_d)./(el_d-el_c);
        r_TVD((el_d-el_c)==0)=1;
        
        % Define Flux Limiter function
        %VANLEER
        phi = (r_TVD + abs(r_TVD))./(1 + abs(r_TVD));
        
        l_TVD=(el_d-el_c)./(el_c-el_up);
        l_TVD((el_c-el_up)==0)=1;
        
        % Define Flux Limiter function
        %VANLEER
        phi_l = (l_TVD + abs(l_TVD))./(1 + abs(l_TVD));        
        
        % Compute fluxes for TVD
        F_rl = a_p.*el_c + a_m.*el_d;
        F_rh = (1/2)*a_vect.*(el_c+el_d) - (1/2)*(a_vect.^2).*...
            (dte./dx_centered).*(el_d-el_c);
        F_ll = a_p.*el_up + a_m.*el_c;
        F_lh= (1/2)*a_vect.*(el_up+el_c) - (1/2)*(a_vect.^2).*...
            (dte./dx_centered).*(el_c-el_up);
             
        F_right = F_rl + phi.*(F_rh-F_rl);
        F_left = F_ll+ phi_l.*( F_lh- F_ll);
        TVD_next= Z(i)-(dte*(F_right-F_left)./dx_centered);
        TVD_next(TVD_next<0)=0;
        
        % In case of nan value, replace by explicit solution        
        exp_next=Z(i)-dte*A(i).*((max(Z(i)-Z(k),0))./dx_ik);        
        TVD_next(isnan(TVD_next))=exp_next(isnan(TVD_next));        
        Z(i) = TVD_next;        
        Z=Z+dte.*upl;
      
end
end
