function [Lx]=derive_L_TVD_x(u_g,a_x,a_p_x,a_m_x,dt,dx,FL,ind)

TVD_center=u_g(ind.inner);
    
    TVD_r_x= u_g(ind.right);
    TVD_r_x2=u_g(ind.right2);
    TVD_l_x= u_g(ind.left);
    
   
    r_TVD_x=zeros(size(a_x));
    r_TVD_x(a_x>0)=(TVD_center(a_x>0)-TVD_l_x(a_x>0))./(TVD_r_x(a_x>0)-TVD_center(a_x>0));
    r_TVD_x(a_x<0)=(TVD_r_x2(a_x<0)-TVD_r_x(a_x<0))./(TVD_r_x(a_x<0)-TVD_center(a_x<0));
    r_TVD_x(diff(TVD_center,[],2)==0)=1;
    
  
    %Define Flux Limiter function
   
   
    phi_x =FL(abs(r_TVD_x));
   
    
    % Compute fluxes for TVD
    F_rl_x = a_p_x(:).*TVD_center + a_m_x(:).*TVD_r_x;
    F_rh_x = (1/2)*a_x(:).*(TVD_center+TVD_r_x) - (1/2)*(a_x(:).^2).*(dt/dx).*(TVD_r_x-TVD_center);
    F_ll_x = a_p_x(:).*TVD_l_x + a_m_x(:).*TVD_center;
    F_lh_x= (1/2)*a_x(:).*(TVD_l_x+TVD_center) - (1/2)*(a_x(:).^2).*(dt/dx).*(TVD_center-TVD_l_x);
    
       % Compute next time step
    phi_x_l=[phi_x(:,1) phi_x(:,1:end-1)];
    F_right_x = F_rl_x + phi_x(:).*(F_rh_x-F_rl_x);
    F_left_x = F_ll_x+ phi_x_l(:).*( F_lh_x- F_ll_x);
    
    Lx=(F_right_x- F_left_x)/dx;
    
    
end