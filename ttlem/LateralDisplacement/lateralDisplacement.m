function [DEM]=lateralDisplacement(DEM,v_x,v_y,time, numerics,BC)
% function to account for horizontal shortening specific to the tectonic
% configuration. Dispalcement is calcualted using 2D numerical advection
% schemes
%
% Syntax
%
%     G = gradient8(DEM,rho_ratio,volume)
%
% Description
%
%     Tectonic shortening is calulcated by imposing an advection field
%     specified in the x and y direction. The phylosophy for this follows
%     from Willet et al.  2001,
%
% Input
%     DEM       digital elevation model (class: GRIDobj)
%     v_x       horizontal shortening in x direction
%     v_x       horizontal shortening in y direction
%     time      time during which displacements take place
%     numerics  properties of numerical methods applied in the model
%               -   numerics.numMeth: 'Upwind_FD' or 'Upwind_TVD'
%               -   numerics.cfl :   courant freidrich lachs number
%     BC        Boundary conditions

% Output
%     DEM       digital elevation model (class: GRIDobj)
%
% Example
%
% See also:
%
% Authors: Benjamin Campforts (benjamin.campforts[at]ees.kuleuven.be)
%          Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% 
% Date: 2 April, 2016

dx=DEM.cellsize;
dy=DEM.cellsize;
u= DEM.Z;
dt=min(dx/max(abs(v_x(:))),dy/max(abs(v_y(:))));
dt=numerics.cfl*dt;
if dt>time
    dt=time;
end
vp_x=max(v_x,0);am_x=min(v_x,0);ap_y=max(v_y,0);am_y=min(v_y,0);

%% Initalisation of main matrix with gohst cells, set boundaries
switch numerics.numMeth
    case 'Upwind_FD'
        BC.nbGhost=1;
    case 'Upwind_TVD'
        BC.nbGhost=2;
end

if BC.nbGhost==1
    u_g=(zeros(size(u,1)+2,size(u,2)+2));
    u_g(2:end-1,2:end-1)=u;
elseif BC.nbGhost==2
    u_g=(zeros(size(u,1)+4,size(u,2)+4));
    u_g(3:end-2,3:end-2)=u;
end
ind=getIndices(u_g,BC);
BC.BC_indices=getBCIndices(BC.nbGhost,u_g);
u_g=setBC(u_g,BC);

%% Run: lateral displacement
for i=1:time/dt
%     disp([num2str(round(i*dt/time*100)) ' % completed'])    
    switch numerics.numMeth
        case 'Upwind_FD'
            %Directional spatial splitting
            Lx=derive_L_FD_x(u_g,vp_x,am_x,ind,dx);
            u_g(ind.inner)=u_g(ind.inner)-dt*Lx;
            u_g=setBC(u_g,BC);
            Ly=derive_L_FD_y(u_g,ap_y,am_y,ind,dy);
            u_g(ind.inner)=u_g(ind.inner)-dt*Ly;
        case 'Upwind_TVD'
            %Directional spatial splitting
            FL = fluxLimiter('superbee');
            Lx=derive_L_TVD_x(u_g,v_x,vp_x,am_x,dt,dx,FL,ind);
            u_g(ind.inner)= u_g(ind.inner)-dt*Lx;
            u_g=setBC(u_g,BC);
            Ly=derive_L_TVD_y(u_g,v_y,ap_y,am_y,dt,dy,FL,ind);
            u_g(ind.inner)= u_g(ind.inner)-dt*Ly;
    end
    u_g=setBC(u_g,BC);
end
DEM.Z=reshape(u_g(ind.inner),DEM.size);
end