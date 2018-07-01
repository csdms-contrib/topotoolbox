function u=Diff_1D(u,D,dx,t_end,BC)
%% @Todo structure BC
%Forward explicit diffusion
% u initial 
% D Diffusion coefficiient in m²/y
% dx spatial resolution
% t time
dt= .9*dx^2/(2*D);
niter=ceil(t_end/dt);
dt=t_end/niter;

[nx]=numel(u);
index=2:nx-1;


for i=1:niter
%     qx=-D*diff(u,1,2)/dx;
%     u(1,2:end-1)=u(1,2:end-1)+(-diff(qx,1,2)/dx)*dt;
    u(index)=u(index)+D*dt/dx^2*(u(index+1)-2*u(index)+u(index-1));
    if isstruct(BC)
        u(1)=BC.start;
        u(end)=BC.end;
    end
end
