function s = sinuosity(S,seglength)

%SINUOSITY sinuosity coefficient 
%
% Syntax
%
%     s = sinuosity(S,seglength)
%
% Description
%
%     The sinuosity is the actual river length divided by the shortest path
%     length. This function calculates the sinuosity for a streamnetwork S
%     along segments of the network with length seglength. seglength is the
%     distance in map units and should be large enough (e.g. 10000 m)
%     depending on the spatial scale of the study.
%
% Input arguments
%
%     S         STREAMobj
%     seglength segment length
%
% Output arguments
%
%     s         node-attribute list with sinuosity
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM);
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S);
%     S = trunk(S);
%     s = sinuosity(S,5000);
%     imageschs(DEM,[],'colormap',[1 1 1],'colorbar',false)
%     hold on
%     plotc(S,s);
%     h = colorbar;
%     h.Label.String = 'Sinuosity';
%
% See also: STREAMobj, STREAMobj/labelreach
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 21. March, 2018


if seglength < (S.cellsize*3)
    error('Choose a larger value of seglength');
end

label = labelreach(S,'seglength',seglength);
d = S.distance;
x = S.x;
y = S.y;
ix = (1:numel(x))';

s = accumarray(label,ix,[max(label) 1],@(ix) sinuosub(ix));
s = s(label);


function sval = sinuosub(ix)

dd = d(ix);
xx = x(ix);
yy = y(ix);

% Euclidean distance
[~,ixmin] = min(dd);
[~,ixmax] = max(dd);
elength = hypot(xx(ixmin)-xx(ixmax),yy(ixmin)-yy(ixmax));

% curvilinear length
clength = dd(ixmax)-dd(ixmin);

% sinuosity
sval  = clength/elength; 

end
end
