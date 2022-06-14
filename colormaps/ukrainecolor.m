function clr = ukrainecolor(n)
%UKRAINECOLOR  Colormap with colors of the Ukrainian flag
%
% Syntax
%
%     y = ukrainecolor(n)
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 28. February, 2022

if nargin == 0
    n = 255;
end

flag = [255   213     0;
        0    91   187];

clr  = interp1([0 255]',double(flag)/255,linspace(0,255,n));



