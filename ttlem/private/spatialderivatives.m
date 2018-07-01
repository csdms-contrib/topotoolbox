function [zx,zy,zxx,zyy,zxy] = spatialderivatives(Z,dx,dy)
% Z needs padding by a 1-pixel wide border

% first derivative kernel
fdk = [0 0 0; -1 0 1; 0 0 0];
zx   = conv2(Z,fdk,'valid')/(2*dx);
zy   = conv2(Z,fdk','valid')/(2*dy);

% second derivative kernel
sdk = [0 0 0; 1 -2 1; 0 0 0];
zxx  = conv2(Z,sdk,'valid')/(dx^2);
zyy  = conv2(Z,sdk','valid')/(dy^2);

zxy  = conv2(Z,[ 1 0 -1; ...
                 0 0  0; ...
                -1 0  1],'valid')/(4*dx*dy); 
end