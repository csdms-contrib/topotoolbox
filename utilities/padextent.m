function padextent(padval,ax)

%PADEXTENT pad current axis extent 
%
% Syntax
%
%     padextent(padval)
%     padextent([padvalx padvaly])
%     padextent([pvleft pvright pvbottom pvtop])
%     
% Description
%
%     padval pads the axes by a specified distance padval. 
%
% Input parameters
%
%     padval     scalar. Will increase the size of the axes in all
%                directions by the specified distance. Negative 
%                values crop the axes.
%     padvalx, padvaly allow you to define different pad distances
%                in x and y direction.
%     pleft, ... let you define values in all directions individually
%     
%
% Example
% 
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     S   = STREAMobj(FD,'minarea',1000);
%     S   = modify(S,'distance',20000);
%     S   = klargestconncomps(S);
%     imageschs(DEM);
%     hold on
%     plot(S,'w')
%     hold off
%     setextent(S)
%     padextent(2000)
%
% See also: IMAGESCHS, GETEXTENT, SETEXTENT
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 23. January, 2020

if nargin == 1
    ax = gca;
end

ext = getextent(ax);

if isscalar(padval)
    ext{1} = ext{1} + [-1 1]*padval;
    ext{2} = ext{2} + [-1 1]*padval;
    setextent(ext,ax);
elseif isnumeric(padval) && numel(padval) == 2
    ext{1} = ext{1} + [-1 1]*padval(1);
    ext{2} = ext{2} + [-1 1]*padval(2);
elseif isnumeric(padval) && numel(padval) == 4
    padval = padval(:)';
    ext{1} = ext{1} + [-1 1].*padval([1 2]);
    ext{2} = ext{2} + [-1 1].*padval([3 4]);
end
setextent(ext,ax);