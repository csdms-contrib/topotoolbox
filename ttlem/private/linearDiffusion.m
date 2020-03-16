function DEM=linearDiffusion(DEM,EYE,p,dx2,C,nrc,L,dt)
% Calcualte linear diffusion using a Crank-Nicolson method
%
% Syntax
%
%       DEM_D=linearDiffusion(DEM,EYE,p,dx2,C,nrc,L,dt)
%
% Description
%
%       Calcualte linear diffusion using a Crank-Nicholson method
%
% Input
%
%       DEM       DEM (digital elevation model) (GRIDobj)
%       EYE       Identity Matrix
%       p         structure array with parameter definitions (see ttlemset)
%       dx2       dx^2
%       C         location of river cells
%       nrc       number of cells in H1.Z
%       L         laplacian
%       dt        time step 
% 
% Output
%
%       DEM       DEM, updated for diffusion (GRIDobj)
%
% Example

if p.AreaThresh > 0;
    D  = EYE + p.D*dt/(2*dx2)*spdiags(~C.Z(:),0,nrc,nrc)*L;
else
    D  = EYE + p.D*dt/(2*dx2)*L;
end
% and solve linear system of equations
% DEM_D = D\DEM.Z(:);
DEM.Z =reshape( pcg_quiet(D,DEM.Z(:),p.DiffTol),DEM.size(1),DEM.size(2));

function [x,flag,relres,iter,resvec] = pcg_quiet(varargin)
[x,flag,relres,iter,resvec] = pcg(varargin{:});
end
end