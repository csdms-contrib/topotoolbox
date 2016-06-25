function [effe] = nonlinear_function(h,delta_x,delta_y,S_c,enne,k,max_nlc)

% NONLINEAR_FUNCTION This function evaluate the nonlinear function
%
% f(h) = k * ( 1 + ( max_nlc -1 ) * ( |grad(h)| / S_c )^n )
%
% k is the diffusivity constant, S_c is the tangent of the critical slope.
% and max_nlc is the maximum value of the non-linear coefficient that
% multiply k.


[h_x,h_y] = gradient(h,delta_x,delta_y);

grad_h = sqrt( h_x.^2 + h_y.^2 );

ratio = min(1.0,grad_h/S_c);

if isinf(enne)
    
    non_linear_coeff = 1.0;
    
else
    
    if ( enne == 0 )
        
        non_linear_coeff = max_nlc;
        
    else
        
        ratio2 = ( ( max_nlc - 1 ) / max_nlc ) ^ ( 1 / enne ) * ratio;
        non_linear_coeff = 1.0 ./ ( 1.0 - ratio2.^enne );
        
    end
    
end


%% wind effect
max_wind_effect = 0.0;

filter_size = 2.0;
h_min = min(min(h));
h_max = max(max(h));


if ( max_wind_effect > 0.0 )
    
    nx = size(h,1);
    ny = size(h,2);
    h_mod = zeros(nx,ny);
    wind_exp = zeros(nx,ny);
    
    angle_rad = 38/180*pi;
    
    for i=1:nx,
        
        h_prof = h(i,:);
        x_prof = (0:1:ny-1)*delta_x;
        
        h_mod(i,:) = h_prof + ( x_prof-x_prof(1) ) * tan(angle_rad);
        
        
        h_wind = 0;
        
        for j = 1:ny,
            
            if  ( h_mod(i,j) >= h_wind )
                
                h_wind = h_mod(i,j);
                wind_exp(i,j) = 1;
                
            else
                
                wind_exp(i,j) = 0;
                
            end
            
        end
        
    end
    
    [h_x,~] = gradient(h_mod,delta_x,delta_y);
    
    %    h_x = abs(min(0.d0,h_x));
    h_x = max(0.01,h_x);
    new_wind = cos(atan(h_x)).^2 .* wind_exp;
    
    
    
    for ij = 1:filter_size,
        
        fls = 1+2*ij;
        
        wind_temp = filter2(ones(fls,fls)/fls^2,new_wind,'valid');
        
        nx2 = (size(wind_exp,1)-size(wind_temp,1))/2;
        ny2 = (size(wind_exp,2)-size(wind_temp,2))/2;
        
        new_wind(1+nx2:nx-nx2,1+ny2:ny-ny2) = wind_temp;
        
    end
    
    new_wind = (1.0/max(max(new_wind))) * new_wind;
    
    new_wind = new_wind .* ( 1 + ((h-h_min)/(h_max-h_min)).^3 );
    
    k = k .* ( 1.0 + max_wind_effect * new_wind );
    
    %     k = k .* max( 0.25 , new_wind.^max_wind_effect );
    
end

% end wind effect

%% sun effect
max_sun_effect = 0.0;

[h_x,h_y] = gradient(h,delta_x,delta_y);


if ( max_sun_effect > 0.d0 )
    
    
    beta = 40;
    chi = 180;
    
    sun_vec(1) = cosd(beta)*cosd(chi);
    sun_vec(2) = cosd(beta)*sind(chi);
    sun_vec(3) = sind(beta);
    
    
    incidence = ( h_x*sun_vec(1) + h_y*sun_vec(2) + 1*sun_vec(3)) ./ ( sqrt( h_x.^2 + h_y.^2 + 1 ));
    
    incidence = incidence - min(min(incidence));
    
    incidence = incidence / max(max(incidence));
    incidence = incidence - mean(mean(incidence));
    
    M = max(max(incidence));
    coeff = ( max_sun_effect - 1 ) / ( M + max_sun_effect * ( 1 - M ) );
    
    %    'coeff'
    %    [max(max(coeff)) min(min(coeff))]
    
    
    sun_effect = 1.0 + coeff * incidence;
    
    %    'sun effect'
    %    [max(max(sun_effect)) min(min(sun_effect))]
    
    
    k = k .* sun_effect;
    
    %    'k'
    %    [max(max(k)) min(min(k))]
    
end


effe = k .* non_linear_coeff;


end
