function Z = funerosion_impnlin(p,Z,dt, A, i,k,dx_ik,upl)

% Implicit solution for nonlinear river incision
%
% Syntax
%
%      Z = funerosion_impnlin(n,Z,dt, A, i,k,dx_ik)
%
%
% Description
%
%       Implicit solution for nonlinear river incision, calucalted by
%       solving the Stream Power Law. This scheme uses a newton rhapson
%       iteration and is unconditionally stable.
%
% Input
%
%       n         slope exponent
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
if p.implCFL
    nrtsteps = ceil(dt/(min(dx_ik)/max(A(:))));
    dte = dt/nrtsteps;
    if dte>3000
        dte=3000;
    end
else
    dte = dt;
end

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
        % z_t
        zt      = Z(i(r));
        % z_(t+dt) of downstream neighbor
        ztp1d   = Z(k(r));
        % dx
        dx      = dx_ik(r);
        
        % initial value for finding root
        if ztp1d < zt;
            ztp1    = newtonraphson(zt,ztp1d,dx,tt,p.n);
        else
            ztp1    = zt;
        end
        
        if ~isreal(ztp1) || isnan(ztp1)
            disp('Non real solutions converted to real')
            ztp1=real(ztp1);
        end
        Z(i(r))=ztp1;
    end
end

    function ztp1 = newtonraphson(zt,ztp1d,dx,tt,n)
        
        tempz   = zt;
        tol = inf;
        
        while tol > 1e-3;
            % iteratively approximated value
            ztp1  =  tempz - (tempz-zt + ...
                (tt*dx) * ...
                ((tempz-ztp1d)./dx)^n) / ...
                (1+n*tt*((tempz-ztp1d)./dx)^(n-1));
            tol   = abs(ztp1-tempz);
            tempz = ztp1;
        end
    end
end