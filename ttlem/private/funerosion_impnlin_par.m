function Z = funerosion_impnlin_par(p,Z,dt, A, S,upl)

% Implicit solution for nonlinear river incision (parallel)
%
% Syntax
%
%      Z = funerosion_impnlin_par(n,Z,dt, A, i,k,dx_ik)
%
%
% Description
%
%       Implicit solution for nonlinear river incision, calculated by
%       solving the Stream Power Law. This scheme uses a newton rhapson
%       iteration and is unconditionally stable. The function runs in
%       parallel and requires the parallel processing toolbox.
%
% Input
%
%       n         slope exponent
%       Z         digital elevation model (matrix)
%       dt        time step
%       A         velocity field
%       S         STREAMobj
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

if isa(S,'STREAMobj');
    CS = STREAMobj2cell(S,'o',20);
else
    S  = STREAMobj(S,'minarea',0);
    CS = STREAMobj2cell(S,'o',20);
end
CZ = cellfun(@(S) Z(S.IXgrid),CS,'UniformOutput',false);
CA = cellfun(@(S) A(S.IXgrid),CS,'UniformOutput',false);

parfor iterpar = 1:numel(CS);
    dx_ik = sqrt((CS{iterpar}.x(CS{iterpar}.ix)-CS{iterpar}.x(CS{iterpar}.ixc)).^2 ...
        + ...
        (CS{iterpar}.y(CS{iterpar}.ix)-CS{iterpar}.y(CS{iterpar}.ixc)).^2);
    CZ{iterpar}  = funerosion_impnlin_par_sub(p.n,CZ{iterpar},dt,CA{iterpar},CS{iterpar}.ix,CS{iterpar}.ixc,dx_ik);
    
end
for iter2 = 1:numel(CS);
    Z(CS{iter2}.IXgrid) = CZ{iter2};
end
end


function z = funerosion_impnlin_par_sub(n,z,dt,a,i,k,dx_ik)

if p.implCFL
    nrtsteps = ceil(dt/(min(dx_ik)/max(A(:))));
else
    nrtsteps=1;
end

% explicit time step
dte = dt/nrtsteps;
time=dt;


while time>0
    time=time-dte;
    if time<0
        dte=dte+time;
        time=0;
    end
    Z=Z+dte.*upl;
    for r = numel(i):-1:1;
        
        tt      = a(i(r))*dte/(dx_ik(r));
        % z_t
        zt      = z(i(r));
        % z_(t+dt) of downstream neighbor
        ztp1d   = z(k(r));
        % dx
        dx      = dx_ik(r);
        
        % initial value for finding root
        if ztp1d < zt;
            ztp1    = newtonraphson(zt,ztp1d,dx,tt,n);
        else
            ztp1    = zt;
        end
        
        if ~isreal(ztp1) || isnan(ztp1)
            disp('Non real solutions converted to real')
            ztp1=real(ztp1);
        end
        z(i(r))=ztp1;
    end
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