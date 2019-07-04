function Z = funerosion_ex(p,Z,dt,dx, A, i,k,dx_ik,upl)
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
%       upl       Vertical uplift of river cells provided as a matrix with
%                 the dimensions of DEM.size
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


% number of time steps
nrtsteps = ceil(dt/(p.cfls*dx/max(A(:))));
% explicit time step
dte = dt/nrtsteps;
time=dt;

a_ori=-A(i);
a=a_ori;
while time>0
    
    if p.n~=1
        s1=((max(Z(i)-Z(k),0))./dx_ik);
        s1(s1<1e-4)=0; %Very low slopes are leading to very small timsteps so they are set to 0.
        exp_f=s1.^(p.n-1);
        exp_f(isinf(exp_f))=1;
        a(exp_f~=0)=a_ori(exp_f~=0).*exp_f(exp_f~=0);
        
        % Check timestep
        dt_calc = p.cfls*min(dx_ik)./max(abs(a));
        if dt_calc<dte
            disp(['For stability, dte is set to: ' num2str(dt_calc)]);
            dte=dt_calc;
        end
    end

    
    time=time-dte;
    if time<0
        dte=dte+time;
        time=0;
    end
    Z=Z+dte.*upl;
    %% Vectorised
%     Z(i) = Z(i) + dte*a.*((max(Z(i)-Z(k),0))./dx_ik);% + U.Z(I)*dte;
    
    %% Not vectorised
    for r = 1:numel(i);   
        Z(i(r))=Z(i(r))-dte*A(i(r)).*((max(Z(i(r))-Z(k(r)),0))./dx_ik(r));        
    end
end