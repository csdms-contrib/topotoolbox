function [r, s] = facet_flow(e0, e1, e2, d1, d2)
%facet_flow Facet flow direction
%
%   [r, s] = facetFlow(e0, e1, e2, d1, d2) computes the facet flow direction
%   for an east-northeast pixel facet.  The facet flow direction is the
%   direction of maximum magnitude slope, constrained to lie within or along
%   an edge of the facet.
%
%   e0 is the height of the center pixel.  e1 is the height of the east
%   neighbor.  e2 is the height of the northeast neighbor.  d1 is the
%   horizontal pixel center spacing.  d2 is the vertical pixel center
%   spacing.  d1 and d2 are optional; if omitted, a value of 1.0 is
%   assumed.
%
%                    e2
%                    |
%                    |
%                    | d2
%                    |
%                    |
%        e0 -------- e1
%              d1
%             
%
%   r is the facet flow direction in radians. r ranges from 0 radians,
%   indicating a flow directly from the center pixel to the east neighbor,
%   to atan(d2, d1) radians, indicating a flow directly from the center
%   pixel to the northeast neighbor.
%
%   s is the downward slope of the facet along the facet flow
%   direction. (Positive s indicates downward flow.) 
%
%   e0, e1, and e2 can be arrays with the same size, in which case r and s
%   are also arrays with the same size.
%
%   Reference: Tarboton, "A new method for the determination of flow
%   directions and upslope areas in grid digital elevation models," Water
%   Resources Research, vol. 33, no. 2, pages 309-319, February 1997.
%
%   Examples
%   --------
%
%      % e1 below e0; e2 below e1. Facet flow is pi/4 radians.
%      e0 = 10; e1 = 5; e2 = 0;
%      [r1, s1] = facet_flow(e0, e1, e2)
%
%      % e2 below e0; e1 below e2. Facet flow is 0 radians.
%      e0 = 10; e2 = 5; e1 = 0;
%      [r2, s2] = facet_flow(e0, e1, e2)
%
%      % e1 and e2 both above e0.  Negative slope value indicates no
%      % downhill flow on the facet.  Returned value for facet flow
%      % direction is not meaningful.  
%      e0 = 0; e1 = 10; e2 = 10;
%      [r3, s3] = facet_flow(e0, e1, e2)
%
%   See also dem_flow, pixel_flow.

%   Steven L. Eddins
%   $Revision: 1.3 $  $Date: 2007/08/02 20:59:13 $
%   Copyright 2007 The MathWorks, Inc.


if nargin < 4
    d1 = 1;
    d2 = 1;
end

s1 = (e0 - e1) / d1;           % eqn (1)
s2 = (e1 - e2) / d2;           % eqn (2)

r = atan2(s2, s1);             % eqn (3)
s = hypot(s1, s2);             % eqn (3)

% Constrain flow direction to lie within or along the edges of the facet.
too_far_south = r < 0;
r(too_far_south) = 0;
s(too_far_south) = s1(too_far_south);

diagonal_angle    = atan2(d2, d1);
diagonal_distance = hypot(d1, d2);
too_far_north = r > diagonal_angle;
r(too_far_north) = diagonal_angle;
s(too_far_north) = (e0(too_far_north) - e2(too_far_north)) / ...
    diagonal_distance;
