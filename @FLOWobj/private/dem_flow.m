function [R, S] = dem_flow(E, d1, d2)
%dem_flow Downslope flow direction for a DEM 
%
%   [R, S] = dem_flow(E, d1, d2) computes the flow direction and downslope
%   for all the pixels in a digital elevation model (E).  E is a matrix of
%   elevation values.  d1 and d2 are the horizontal and vertical pixel
%   center spacing.  d1 and d2 are optional; if omitted, a value of 1.0 is
%   assumed.
%
%   R, a matrix the same size as E, contains the pixel flow direction, in
%   radians, for each pixel of E.  Pixel flow direction is measured counter
%   clockwise from the east-pointing horizontal axis.  R is NaN for each pixel
%   of E that has no downhill neighbors.
%
%   S, a matrix the same size as E, contains the downward slope (along the
%   pixel flow direction) for each pixel of E.  Negative values indicate that
%   the corresponding pixel of E has no downhill neighbors.
%
%   Reference: Tarboton, "A new method for the determination of flow
%   directions and upslope areas in grid digital elevation models," Water
%   Resources Research, vol. 33, no. 2, pages 309-319, February 1997. 
%
%   Example
%   -------
%
%       s = load('milford_ma_dem');
%       E = s.Zc(146:247, 146:210);
%       [R, S] = dem_flow(E);
%       vis_dem_flow(E, R, S);
%
%   See also facet_flow, pixel_flow, upslope_area, vis_dem_flow.

%   Steven L. Eddins
%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.4 $  $Date: 2007/08/02 20:59:13 $

if nargin < 3
    d1 = 1;
    d2 = 1;
end

% Pad E so we can calculate pixel flow along the DEM boundaries.
Ep = padarray(E, [1 1], 'replicate', 'both');
[M,N] = size(Ep);

[i, j] = ndgrid(2:M-1, 2:N-1);
[R, S] = pixel_flow(Ep, i, j, d1, d2);
