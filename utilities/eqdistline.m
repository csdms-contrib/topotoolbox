function [newx,newy] = eqdistline(x,y,step)

% Function for getting a new line with equi-distant (step) vertices from 
% irregular line. Assumes that x,y coordinates are in same units as step.
% Code as explained in:
% http://blogs.mathworks.com/steve/2012/07/06/walking-along-a-path/

dx = diff(x);
dy = diff(y);
dxy = hypot(dx, dy); % square root of sum of squares
cumdxy = [0; cumsum(dxy,1)];

% get rid of duplicates
repeatedPtsMask = diff(cumdxy)==0;
if any(repeatedPtsMask)
    x(repeatedPtsMask) = [];
    y(repeatedPtsMask) = [];
    cumdxy(repeatedPtsMask) = [];
end

npoints = 1+floor(cumdxy(end)/step);
dist_steps = linspace(0, cumdxy(end), npoints);

newx = interp1(cumdxy, x, dist_steps);
newy = interp1(cumdxy, y, dist_steps);