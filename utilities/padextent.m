function extout = padextent(padval,varargin)

%PADEXTENT pad current axis extent 
%
% Syntax
%
%     padextent(padval)
%     padextent([padvalx padvaly])
%     padextent([pvleft pvright pvbottom pvtop])
%     padextent(...,ax)
%     padextent(...,pn,pv)
%     extent = padextent(...)
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
%     ax         handle to axes (by default gca)
%
%     Parameter name/value pairs
%
%     'unit'     {'map'} or 'percent'
%     
% Example 1
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
% Example 2: Zoom in
%
%     imageschs(DEM,[],'colorbar',false,'ticklabels','none')
%     for r = 1:100; padextent(-1,'unit','perc'); pause(.01); end
%
% See also: IMAGESCHS, GETEXTENT, SETEXTENT
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 23. January, 2020

p = inputParser;
addRequired(p,'padval');
addOptional(p,'ax',gca)
addParameter(p,'unit','map')
parse(p,padval,varargin{:});

ax  = p.Results.ax;
ext = getextent(ax);

unit = validatestring(p.Results.unit,{'map' 'percent'});

% padextent will use map units, hence any conversion is done here
switch unit
    case 'percent'
        padval = padval/100;
        if isscalar(padval)
            padvalmapx = padval * abs(diff(ext{1}));
            padvalmapy = padval * abs(diff(ext{2}));
            padval = [padvalmapx padvalmapy];
        elseif isnumeric(padval) && numel(padval) == 2
            padvalmapx = padval(1) * abs(diff(ext{1}));
            padvalmapy = padval(2) * abs(diff(ext{2}));
            padval = [padvalmapx padvalmapy];
        elseif isnumeric(padval) && numel(padval) == 4
            padvalmapx = padval([1 2]) * abs(diff(ext{1}));
            padvalmapy = padval([3 4]) * abs(diff(ext{2}));
            padval = [padvalmapx padvalmapy];
        end
end



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

if nargout == 1
    extout = ext;
end