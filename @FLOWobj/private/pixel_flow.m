function [R, S] = pixel_flow(E, i, j, d1, d2)
%pixel_flow Downslope flow direction for DEM pixels
%
%   [R, S] = pixel_flow(E, i, j, d1, d2) computes the flow direction and
%   slope for specified pixels (i and j) of a digital elevation model (E).  E
%   is a matrix of elevation values.  i and j, which must be the same size,
%   are the row and column coordinates of the specified pixels.  d1 and d2
%   are the horizontal and vertical pixel center spacing.  d1 and d2 are
%   optional; if omitted, a value of 1.0 is assumed.
%
%   The specified pixels cannot be on the border of E.  In other words, i
%   must be greater than 1 and less than size(E, 1), and j must be greater
%   than 1 and less than size(E, 2).
%
%   R, which is the same size as i and j, contains the pixel flow direction
%   in radians.  Pixel flow direction is measured counter clockwise from the
%   east-pointing horizontal axis.  R is NaN for each pixel that has no
%   downhill neighbors. 
%
%   S, which is the same size as i and j, is the downward slope along the
%   pixel flow direction.  Negative values indicate that the corresponding
%   pixels have no downhill neighbors.
%
%   Note: Connected groups of NaN pixels touching the border are treated as
%   being at a higher elevation than all the other pixels in E.
%
%   Reference: Tarboton, "A new method for the determination of flow
%   directions and upslope areas in grid digital elevation models," Water
%   Resources Research, vol. 33, no. 2, pages 309-319, February 1997. 
%
%   Examples
%   --------
%
%      % Flow from the center pixel goes to the right, so the pixel flow
%      % direction is 0 radians.
%      E1 = [2 1 0; 2 1 0; 2 1 0]
%      R1 = pixel_flow(E1, 2, 2)
%
%      % Flow from the center pixel goes to the upper right, so the pixel
%      % flow direction is pi/4 radians.
%      E2 = [2 1 0; 3 2 1; 4 3 2]
%      R2 = pixel_flow(E2, 2, 2)
%
%      % The center pixel has no downhill neighbors, so the pixel flow
%      % direction is NaN.
%      E3 = [2 2 2; 2 1 2; 2 2 2]
%      R3 = pixel_flow(E3, 2, 2)
%
%   See also dem_flow, facet_flow.

%   Steven L. Eddins
%   $Revision: 1.5 $  $Date: 2007/10/02 15:50:06 $
%   Copyright 2007 The MathWorks, Inc.


if nargin < 4
    d1 = 1;
    d2 = 1;
end

[M, N] = size(E);

% Preprocess NaNs connected to the border.  Make them higher than the
% highest non-NaN value.
highest_value = max(E(:));
bump = min(1, eps(highest_value));
E(border_nans(E)) = highest_value + bump;

% Compute linear indices at desired locations.
e0_idx = (j - 1)*M + i;

% Table 1, page 311
% Row and column offsets corresponding to e1 and e2 for each
% table entry:
e1_row_offsets = [0 -1 -1  0  0  1  1  0];
e1_col_offsets = [1  0  0 -1 -1  0  0  1];

e2_row_offsets = [-1 -1 -1 -1  1  1  1  1];
e2_col_offsets = [ 1  1 -1 -1 -1 -1  1  1];


% Linear e1 and e2 offsets.
e1_linear_offsets = e1_col_offsets*M + e1_row_offsets;
e2_linear_offsets = e2_col_offsets*M + e2_row_offsets;


% Initialize R and S values based on the first facet.
E0 = E(e0_idx);
E1 = E(e0_idx + e1_linear_offsets(1));
E2 = E(e0_idx + e2_linear_offsets(1));


[R, S] = facet_flow(E0, E1, E2, d1, d2);

% Multipliers ac and af are used to convert a flow angle from the canonical
% east-northeast facet to any of the eight facets.
ac = [0  1  1  2  2  3  3  4];
af = [1 -1  1 -1  1 -1  1 -1];

positive_S = S > 0;
R(positive_S) = (af(1) * R(positive_S)) + (ac(1) * pi / 2);  % Equation (6)
R(~positive_S) = NaN;

for k = 2:8
    % Compute Rk and Sk corresponding to the k-th facet. Where Sk is positive
    % and greater than S, replace S and recompute R based on Rk.
    E1 = E(e0_idx + e1_linear_offsets(k));
    E2 = E(e0_idx + e2_linear_offsets(k));
    [Rk, Sk] = facet_flow(E0, E1, E2, d1, d2);
    
    new_R = (Sk > S) & (Sk > 0);
    S = max(S, Sk);
    
    R(new_R) = (af(k) * Rk(new_R)) + (ac(k) * pi / 2);
end
