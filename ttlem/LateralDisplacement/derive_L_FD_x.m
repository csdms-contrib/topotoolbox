function Lx=derive_L_FD_x(u_g,ap_x,am_x,ind,dx)
ux_m=(u_g(ind.inner)-u_g(ind.left))/dx;
ux_p=(u_g(ind.right)-u_g(ind.inner))/dx;
Lx=(ap_x(:).*ux_m+am_x(:).*ux_p);
end