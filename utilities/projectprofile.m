function [d,z]=projectprofile(x,y,z)

% 
%
% [d,znew]=projectprofile(x,y,z)
% _______________________________________________________
%
% When surveying a cross profile using DGPS, projectprofile 
% helps you to visualize the profile by projecting each 
% measured point on the straight line between the starting 
% and the end point of the profile.
%
% x, y and z must be column vectors of same length.
% _______________________________________________________
% WS
%


% starting point
x = x(:);
y = y(:);
z = z(:);

if x(1) > x(end);
    x = flipud(x);
end
sp = [x(1); y(1)];

% subtract starting point from other locations
x = x-sp(1);
y = y-sp(2);

% end point
a = [x(end); y(end)];

% projection of [x; y] on vector a
beta = (a' * [x y]')./(a'*a);
d = beta .* sqrt(a'*a);

% sort the distance vector and z ascending
% [d,IX1] = sort(d);
% znew = z(IX1);

% minimum point is zero
% [zmin,IX2] = min(znew);
% d = d-d(IX2(1));
% d = d';

