function Ly=derive_L_FD_y(u_g,ap_y,am_y,ind,dy)
uy_m=(u_g(ind.inner)-u_g(ind.top))/dy;
uy_p=(u_g(ind.bottom)-u_g(ind.inner))/dy;
Ly=(ap_y(:).*uy_m+am_y(:).*uy_p);
end