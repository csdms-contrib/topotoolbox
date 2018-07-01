function DEM = shortening(DEM,v_x,v_y,dt)

% function to account for horizontal shortening specific to the tectonic configuration 
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
%
%     DEM       digital elevation model (class: GRIDobj)
%     v_x       horizontal shortening in x direction
%     v_x       horizontal shortening in y direction
%     dt        timestep
%
% Output
%
%     DEM       digital elevation model (class: GRIDobj)
%                  
% Example
% 
% See also: 
% 
% Authors: Benjamin Campforts (benjamin.campforts[at]ees.kuleuven.be)
%          Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
%
% Date: 25 March, 2016
%Check Courant
dt_cal=min(abs((0.95./((v_x(:)+v_y(:))./DEM.cellsize))));

nrtsteps = ceil(dt/dt_cal);
% explicit time step
dte = dt/nrtsteps;
time=dt;

while time>0;
    time=time-dte;
    if time<0
        dte=dte+time;
        time=0;
    end
    change=zeros(size(DEM.Z));
    
    change(:,1:end-1) =change(:,1:end-1) - (v_x(:,2:end)<0).*v_x(:,2:end).*diff(DEM.Z,1,2)/DEM.cellsize;     % advection
    change(:  ,2:end) =change(:  ,2:end) - (v_x(:,1:end-1)>0).*v_x(:,1:end-1).*diff(DEM.Z,1,2)/DEM.cellsize;
    change(1:end-1,:) =change(1:end-1,:) - (v_y(2:end,:)<0).*v_y(2:end,:).*diff(DEM.Z,1,1)/DEM.cellsize;
    change(2:end,:) = change(2:end,:)- (v_y(1:end-1,:)>0).*v_y(1:end-1,:).*diff(DEM.Z,1,1)/DEM.cellsize;
    
    DEM.Z=DEM.Z+change*dte;    
    
end
%%
