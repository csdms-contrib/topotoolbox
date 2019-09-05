function c = cumsum(S,c,direction)

%CUMSUM cumulative sum on stream network
%
% Syntax
%
%     c = cumsum(S,a)
%     c = cumsum(S,a,direction)
%
% Description
%
%     cumsum calculates the cumulative sum along a stream network. By
%     default, cumsum accumulates the values in downstream direction.
%
% Input arguments
%
%     S          STREAMobj
%     a          node-attribute list
%     direction  'downstream' (default) or 'upstream'
%
% Output
%
%     c    node attribute list (double)
%
% Example
%
%     % Two step calculation of upslope area
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     % flow direction
%     FD  = FLOWobj(DEM,'preprocess','carve','mex',true);
%     % stream network
%     S   = STREAMobj(FD,'minarea',500);
%     a1  = hillslopearea(S,FD);
%     a2  = cumsum(S,a1);
%
% See also: STREAMobj/cumtrapz
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 19. August 2019

narginchk(2,3)

if isa(c,'GRIDobj')
    c = getnal(S,c);
else
    if ~isnal(S,c)
        error('The second input argument must be a node attribute list.')
    end
end

c = double(c);

if nargin <= 2
    direction = 'downstream';
else 
    direction = validatestring(direction,{'upstream','downstream'},'STREAMobj/cumsum','direction',3);
end

% get sorted edge list 
ix  = S.ix;
ixc = S.ixc;

% traverse stream network 
switch direction
    case 'upstream'
        for r = numel(S.ix):-1:1
            c(ix(r)) = c(ixc(r)) + c(ix(r));
        end
    otherwise
        for r = 1:numel(S.ix)
            c(ixc(r)) = c(ix(r)) + c(ixc(r));
        end
end
