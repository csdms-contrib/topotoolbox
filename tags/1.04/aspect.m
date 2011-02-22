function ASP = aspect(X,Y,dem)

% derive aspect from a digital elevation model
%
% Syntax
%
%     ASP = aspect(X,Y,dem)
%
% Description
%
%     aspect returns the slope exposition of each cell in a digital
%     elevation model in degrees. In contrast to the second output of
%     gradient8 which returns the steepest slope direction, aspect 
%     returns the angle as calculated by surfnorm.
%
% Input
%
%     X,Y       coordinate matrices as generated by meshgrid
%     dem       digital elevation model
%
% Output
% 
%     ASP       Aspect in degrees (clockwise from north)
%
% Example
%
%     load exampleDEM
%     ASP = aspect(X,Y,dem);
%     surf(X,Y,dem,ASP)
%
%
% See also: gradient8
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 24. February, 2009


if Y(1)>Y(2)
    flagflip = false;
else
    Y = flipud(Y);
    dem = flipud(dem);
    flagflip = true;
end


[Nx,Ny] = surfnorm(X,Y,dem);

ASP = cart2pol(Ny,Nx);
ASP = ASP/pi * 180 + 180;

if flagflip
    ASP = flipud(ASP);
end