function [Ly]=derive_L_TVD_y(u_g,a_y,a_p_y,a_m_y,dt,dy,FL,ind)

TVD_center=u_g(ind.inner);

TVD_r_y= u_g(ind.bottom);
TVD_r_y2=u_g(ind.bottom2);
TVD_l_y= u_g(ind.top);

r_TVD_y=zeros(size(a_y));
r_TVD_y(a_y>0)=(TVD_center(a_y>0)-TVD_l_y(a_y>0))./(TVD_r_y(a_y>0)-TVD_center(a_y>0));
r_TVD_y(a_y<0)=(TVD_r_y2(a_y<0)-TVD_r_y(a_y<0))./(TVD_r_y(a_y<0)-TVD_center(a_y<0));
r_TVD_y(abs(TVD_r_y-TVD_center)<2*eps)=0;
phi_y =FL(abs(r_TVD_y));

% Compute fluxes for TVD
F_rl_y = a_p_y(:).*TVD_center + a_m_y(:).*TVD_r_y;
F_rh_y = (1/2)*a_y(:).*(TVD_center+TVD_r_y(:)) - (1/2)*(a_y(:).^2).*(dt/dy).*(TVD_r_y(:)-TVD_center);
F_ll_y = a_p_y(:).*TVD_l_y + a_m_y(:).*TVD_center;
F_lh_y= (1/2)*a_y(:).*(TVD_l_y+TVD_center) - (1/2)*(a_y(:).^2).*(dt/dy).*(TVD_center-TVD_l_y);

% Compute next time step
phi_y_l=[phi_y(1,:); phi_y(1:end-1,:)];
%     phi_y_l=[phi_x(:,1) phi_x(:,1:end-1)];
F_right_y = F_rl_y +phi_y(:).*(F_rh_y-F_rl_y);
F_left_y = F_ll_y+ phi_y_l(:).*( F_lh_y- F_ll_y);
Ly=(F_right_y- F_left_y)/dy;

end