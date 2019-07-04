function u=Diff_2D(u,D,dx,dy,t_end,BC)
%Forward explicit diffusion
% u initial
% D Diffusion coefficiient in m²/y
% dx spatial resolution in x
% dx spatial resolution in y
% t time
%@ Example
% BCDiff.nbGost=1;
% BCDiff.type='Dirichlet_Matrix';
% BCDiff.dir_value=H1.Z;
% BCDiff=getBCIndices(BCDiff,H1.Z);%
% Z=Diff_2D(H1.Z,p.D,H1.cellsize,H1.cellsize,dt,BCDiff);
% HD=Z(:);
 %-----------------------------------------------------------
 
 
if numel(D)==1
    D_x=D;
    D_y=D;
else
    D_x=D(1);
    D_y=D(2);
end

dt=.9*dx^2*dy^2/(2*D_x*dx^2+2*D_y*dy^2);
niter=ceil(t_end/dt);
dt=t_end/niter;


[nx,ny]=size(u);
ind_i=2:nx-1;
ind_j=2:ny-1;
for i=1:niter
    u(ind_i,ind_j)=u(ind_i,ind_j)+(D_y*dt*(u(ind_i+1,ind_j)-2*u(ind_i,ind_j)+u(ind_i-1,ind_j))/(dy^2))+(D_x*dt*(u(ind_i,ind_j+1)-2*u(ind_i,ind_j)+u(ind_i,ind_j-1))/(dx^2));
end


switch BC.type
    case 'Dirichlet'
        u([BC.BC_indices.topRow BC.BC_indices.botRow BC.BC_indices.leftRow BC.BC_indices.rightRow])=BC.dir_value;
    case 'Dirichlet_Matrix'
        u([BC.BC_indices.topRow BC.BC_indices.botRow BC.BC_indices.leftRow BC.BC_indices.rightRow])=BC.dir_value([BC.BC_indices.topRow BC.BC_indices.botRow BC.BC_indices.leftRow BC.BC_indices.rightRow]);
end
end