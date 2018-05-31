function Z = funerosion_implin(p,Z,dt, A, i,k,dx_ik,upl)

% Implicit solution for river incision
%
% Syntax
%
%      Z = funerosion_implin(Z,dt,A, i,k,dx_ik)
%
%
% Description
%
%       Implicit solution for river incision, calucalted by solving the Stream Power
%       Law. This scheme is unconditionally stable.
%
% Input
%
%       Z         digital elevation model (matrix)
%       dt        time step
%       A         velocity field
%       i         giver, derived from the flow direction FD.ix where ix is an edge attribute and represetns topologically sorted nodes (givers)
%       k         receiver, derived from the flow direction FD.ixc where ixc is an edge attribute and represetns topologically sorted nodes (receivers)
%       dx_ik     horizontal distance
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

if p.implCFL
    nrtsteps = ceil(dt/(min(dx_ik)/max(A(:))));
    dte = dt/nrtsteps;
    if dte>3000
        dte=3000;
    end
else
    dte = dt;
end

% explicit time step

time=dt;


while time>0
    time=time-dte;
    if time<0
        dte=dte+time;
        time=0;
    end
    Z=Z+dte.*upl;    
    
    for r = numel(i):-1:1;
        tt      = A(i(r))*dte/(dx_ik(r));
        Z(i(r)) = (Z(i(r)) + Z(k(r))*tt)./(1+tt);
    end
end
end