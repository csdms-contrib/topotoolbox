function [flux_east,flux_west,flux_north,flux_south] = ...
    numerical_fluxes( h , f , delta_x , delta_y )

% NUMERICAL_FLUXES This function evaluates the fluxes at the boundaries
% of each cell (east,west,south,north) with a finite-difference scheme.
 

nx = size(h,1);
ny = size(h,2);

flux_east = zeros(nx,ny);
flux_west = zeros(nx,ny);
flux_north = zeros(nx,ny);
flux_south = zeros(nx,ny);


flux_west(2:nx,:) = 0.5D0 * ( f(1:nx-1,:) + f(2:nx,:) ) .* ...
    ( h(2:nx,:) - h(1:nx-1,:) ) / delta_x ;

flux_east(1:nx-1,:) = flux_west(2:nx,:);

flux_south(:,2:ny) = 0.5D0 * ( f(:,1:ny-1) + f(:,2:ny) ) .* ...
    ( h(:,2:ny) - h(:,1:ny-1) ) / delta_y ;

flux_north(:,1:ny-1) = flux_south(:,2:ny);

%% Attempt to limit the fluxes

% max_flux(1:nx,1:ny) = h_avl(1:nx,1:ny) * delta_x * delta_y / delta_t;
% 
% flux_west_new = flux_west;
% 
% flux_west_new(2:nx,1:ny) = ( flux_west(2:nx,1:ny) > 0 ) ...
%     .* min( flux_west(2:nx,1:ny) , max_flux(1:nx-1,1:ny) ) ...
%     + ( flux_west(2:nx,1:ny) < 0 ) ...
%     .* max( flux_west(2:nx,1:ny) , - max_flux(2:nx,1:ny) );
% 
% flux_east_new = flux_east;
% 
% flux_east_new(1:nx-1,1:ny) = ( flux_east(1:nx-1,1:ny) > 0 ) ...
%     .* min( flux_east(1:nx-1,1:ny) , max_flux(1:nx-1,1:ny) ) ...
%     + ( flux_east(1:nx-1,1:ny) < 0 ) ...
%     .* max( flux_east(1:nx-1,1:ny) , - max_flux(2:nx,1:ny) );
% 
% flux_south_new = flux_south;
% 
% flux_south_new(1:nx,2:ny) = ( flux_south(1:nx,2:ny) > 0 ) ...
%     .* min( flux_south(1:nx,2:ny) , max_flux(1:nx,1:ny-1) ) ...
%     + ( flux_south(1:nx,2:ny) < 0 ) ...
%     .* max( flux_south(1:nx,2:ny) , - max_flux(1:nx,2:ny) );
% 
% flux_north_new = flux_north;
% 
% flux_north_new(1:nx,1:ny-1) = ( flux_north(1:nx,1:ny-1) > 0 ) ...
%     .* min( flux_north(1:nx,1:ny-1) , max_flux(1:nx,1:ny-1) ) ...
%     + ( flux_north(1:nx,1:ny-1) < 0 ) ...
%     .* max( flux_north(1:nx,1:ny-1) , - max_flux(1:nx,2:ny) );
% 
% flux_west = flux_west_new;
% flux_east = flux_east_new;
% flux_south = flux_south_new;
% flux_north = flux_north_new;

% min f*x such that A*x<=0 , Aeq * x = beq , lb <= x <= ub
% Aeq = [];
% beq = [];
% lb = 0.0;
% ub = 1.0;
% f = 1.0;

% x = linprog(f,A,b,Aeq,beq,lb,ub,x0)

end
