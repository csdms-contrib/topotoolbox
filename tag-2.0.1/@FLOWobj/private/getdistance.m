function d = getdistance(ix,ixc,siz,cellsize)

% distance between each node
%
% Syntax
%
%     d = getdistance(ix,ixc,siz,cellsize)
%
% Description
%
%     calculate the distance between each node which is either
%     sqrt(2)*cellsize or cellsize.
%
% 
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 26. January, 2013

d = ones(size(ix))*cellsize;
diagdistance = norm([cellsize cellsize]);
nrrows = siz(1);
if mod(nrrows,2) == 1
    % odd number of rows
    d(mod(ix,2) == mod(ixc,2)) = diagdistance;
else
    % even number of rows 
    d(mod(ix,2) ~= mod(ixc,2) & ~(imabsdiff(ix,ixc)==1)) = diagdistance; % ~(ix-ixc == 1 | ixc-ix == 1))
end
    