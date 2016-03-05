function [newx,newy,newdist] = interpline(x,y,d,step)

% distribute vertices evenly along line with irregularly spaced vertices
%
% Syntax
%
%     [newx,newy,newdist] = interpline(x,y,d,step)
%
% Input
%
%     x, y       coordinate vectors
%     d          cumulative distance along path defined by x and y
%     step       new distance between vertices
%
%
% Function for getting a new line with equi-distant (step) vertices from
% irregular line. Assumes that x,y coordinates are in same units as step.
% Code as explained in:
% http://blogs.mathworks.com/steve/2012/07/06/walking-along-a-path/
if isrow(x); x = x'; end
if isrow(y); y = y'; end
if isrow(d); d = d'; end
% get rid of duplicates
repeatedPtsMask = diff(d)==0;
if any(repeatedPtsMask)
    x(repeatedPtsMask) = [];
    y(repeatedPtsMask) = [];
    d(repeatedPtsMask) = [];
end
% interpolate
newdist = (d(1):step:d(end))';
newx = interp1(d, x, newdist);
newy = interp1(d, y, newdist);
end %