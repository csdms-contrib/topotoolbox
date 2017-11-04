function T = flow_matrix(E, R, d1, d2)
%flow_matrix System of linear equations representing pixel flow
%
%   T = flow_matrix(E, R) computes a sparse linear system representing flow from
%   pixel to pixel in the DEM represented by the matrix of height values, E.  R
%   is the matrix of pixel flow directions as computed by dem_flow.  T is
%   numel(E)-by-numel(E).  The value T(i,j) is the negative of the fractional
%   flow from pixel j to pixel i, where pixels are numbered columnwise. For
%   example, if E is 15-by-10, then T is 150-by-150, and T(17,18) is the
%   negative of the fractional flow from pixel 18 (row 3, column 2) to pixel 17
%   (row 2, column 2).
%
%   d1 is the horizontal pixel spacing.  d2 is the vertical pixel spacing.  d1
%   and d2 are optional; if omitted, a value of 1.0 is assumed.
%
%   Note: Connected groups of NaN pixels touching the border are treated as
%   having no contribution to flow.
%
%   Reference: Tarboton, "A new method for the determination of flow
%   directions and upslope areas in grid digital elevation models," Water
%   Resources Research, vol. 33, no. 2, pages 309-319, February 1997. 
%
%   Example
%   -------
%
%       E = peaks;
%       R = dem_flow(E);
%       T = flow_matrix(E, R);
%       spy(T)
%
%   See also dem_flow, dependency_map, influence_map, upslope_area.

%   Steven L. Eddins
%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.4 $  $Date: 2008/02/14 15:49:21 $

if nargin < 4
    d1 = 1;
    d2 = 1;
end

[M, N] = size(R);
nrc    = numel(R);

% Label the pixels.
pixel_labels = reshape(1:numel(R), M, N);

% What are the rows and columns of the pixels that have a north neighbor?
have_neighbor.north.rows = 2:M;
have_neighbor.north.cols = 1:N;

% What are the rows and columns of the pixels that have a northeast
% neighbor?
have_neighbor.northeast.rows = 2:M;
have_neighbor.northeast.cols = 1:N-1;

% And so on...
have_neighbor.east.rows = 1:M;
have_neighbor.east.cols = 1:N-1;

have_neighbor.southeast.rows = 1:M-1;
have_neighbor.southeast.cols = 1:N-1;

have_neighbor.south.rows = 1:M-1;
have_neighbor.south.cols = 1:N;

have_neighbor.southwest.rows = 1:M-1;
have_neighbor.southwest.cols = 2:N;

have_neighbor.west.rows = 1:M;
have_neighbor.west.cols = 2:N;

have_neighbor.northwest.rows = 2:M;
have_neighbor.northwest.cols = 2:N;

% What are the linear index offsets corresponding to each of the neighbors?
offset.north     = -1;
offset.northeast = M - 1;
offset.east      = M;
offset.southeast = M + 1;
offset.south     = 1;
offset.southwest = -M + 1;
offset.west      = -M;
offset.northwest = -M - 1;

directions = {'north', 'northeast', 'east', 'southeast', ...
    'south', 'southwest', 'west', 'northwest'};

% Initialize the index and value vectors that will be used to construct the
% sparse array.
ii = [];
jj = [];
vv = [];

for k = 1:numel(directions)
    direction = directions{k};
    
    i = pixel_labels(have_neighbor.(direction).rows, ...
        have_neighbor.(direction).cols);
    j = i + offset.(direction);
    p = -directional_weight(direction, R(j), d1, d2);
    
    non_zero_weight = p ~= 0;
    
    ii = [ii; i(non_zero_weight)]; %#ok<AGROW>
    jj = [jj; j(non_zero_weight)]; %#ok<AGROW>
    vv = [vv; p(non_zero_weight)]; %#ok<AGROW>
end
    
% Filter out weights calculated for border NaN pixels.
border_mask = border_nans(E);
border_labels = pixel_labels(border_mask);
delete_mask = ismember(ii, border_labels) | ...
    ismember(jj, border_labels);
ii(delete_mask) = [];
jj(delete_mask) = [];
vv(delete_mask) = [];

% The weight for the pixel itself is 1.
% ii = [ii; pixel_labels(:)];
% jj = [jj; pixel_labels(:)];
% vv = [vv; ones(numel(R), 1)];

[ip, jp, vp] = plateau_flow_weights(E, R, border_mask);

ii = [ii; ip];
jj = [jj; jp];
vv = [vv; vp];

% Construct the sparse array representing the system of linear pixel flow
% equations.

T = sparse(ii, jj, double(vv),nrc,nrc);

%==========================================================================
function [ip, jp, vp] = plateau_flow_weights(E, R, border_mask)
%
% For the pixels with no assigned flow direction (marked in R with NaN), compute
% how their flow is to be distributed to their neighbors.  The output vectors
% (ip, jp, and vp) are intended to be used in the formation of a sparse linear
% system.  T(i, j) = v indicates that the fraction -v of the flow from pixel j
% should be directed to pixel i.
%
% Pixels with no flow direction that are part of regional minima do not
% distribute any flow to their neighbors.
%
% The algorithm was inspired by the "arrowing" algorithm described in
% F. Meyer, "Skeletons and Perceptual Graphs," Signal Processing 16 (1989)
% 335-363.  See Appendix A.2, Arrowing, on pages 360-361.
%
% Flow allocations are computed iteratively.  In each iteration, find all the
% pixels that have no computed flow, and that have equal-valued neighbors
% with a previously computed flow.  Distribute pixel flow equally among these
% neighbors.  Continue until all pixels have been assigned a flow
% distribution, except for those belonging to regional minima.

S = size(R);

% Replace border NaNs with Inf.  imregionalmin, used below, doesn't like
% NaNs in its input.
E(border_mask) = Inf;

% Initialize list of pixels whose flow is marked NaN. Don't include pixels
% that belong to regional minima, and don't include pixels that are
% border NaNs in E.
[nan_list_r, nan_list_c] = find(isnan(R) & ~imregionalmin(E) & ...
    ~border_mask);

done = false;

num_nans = numel(nan_list_r);
ip = zeros(8*num_nans, 1);
jp = zeros(8*num_nans, 1);
vp = zeros(8*num_nans, 1);
total_count = 0;

while ~done
    done = true;

    % Use this variable to keep track of pixels that won't need to be visited
    % in the next iteration of the while loop.
    delete_from_list = false(numel(nan_list_r), 1);

    % Use rr and cc to track the locations in R that should be set to 0 at
    % the end of this loop iteration.
    rr = zeros(numel(nan_list_r), 1);
    cc = zeros(numel(nan_list_c), 1);
    
    % Use ww and zz to track the set of neighbors that should receive
    % flow from a given pixel.
    ww = zeros(8,1);
    zz = zeros(8,1);
    nan_count = 0;
    for k = 1:numel(nan_list_r)
        r = nan_list_r(k);
        c = nan_list_c(k);

        % Find all neighbors (w,z) of (r,c) for which R(w,z) is not NaN, and
        % for which E(w,z) = E(r,c).  Count how many there are and record
        % their locations.  Pixel flow will be distributed to all such
        % neighbors equally.
        neighbor_count = 0;
        for w = max(1, r-1) : min(S(1), r+1)
            for z = max(1, c-1) : min(S(2), c+1)
                if ~isnan(R(w,z)) && (E(r,c) == E(w,z))
                    neighbor_count = neighbor_count + 1;
                    ww(neighbor_count) = w;
                    zz(neighbor_count) = z;
                end
            end
        end
        
        if neighbor_count > 0
            % We haven't converged.  At least one more execution of the while
            % loop will be required.
            done = false;
            
            nan_count = nan_count + 1;
            rr(nan_count) = r;
            cc(nan_count) = c;

            % Compute indices and values to be added to ip, jp, and vp.  i
            % contains the linear indices of the neigbhors of (r,c).  j
            % contains the linear index of (r,c), copied neighbor_count
            % times.  v is calculated so that each neighbor gets an equal
            % share of the pixel flow from (r,c).
            i = (zz(1:neighbor_count) - 1)*S(1) + ww(1:neighbor_count);
            j = ((c - 1)*S(1) + r) * ones(neighbor_count, 1);
            v = -ones(neighbor_count, 1) / neighbor_count;
            
            slots = total_count + (1:neighbor_count);
            ip(slots) = i;
            jp(slots) = j;
            vp(slots) = v;
            total_count = total_count + neighbor_count;
            
            delete_from_list(k) = true;
        end
    end
    rr = rr(1:nan_count);
    cc = cc(1:nan_count);
    
    if ~done
        nan_list_r(delete_from_list) = [];
        nan_list_c(delete_from_list) = [];
        
        R((cc - 1)*S(1) + rr) = 0;
    end
end

ip = ip(1:total_count);
jp = jp(1:total_count);
vp = vp(1:total_count);
%--------------------------------------------------------------------------

%==========================================================================
function w = directional_weight(direction, R, d1, d2)
%
% Given a direction (from pixel to neighbor), and an array of neighbor flow
% directions, R, compute the corresponding pixel flow weights.
%
% See column 1, page 313 of the Tarboton paper.

% Arrange the angles from a pixel to its neighbor in a counterclockwise
% (increasing angle) sequence.
theta_d = atan2(d2, d1);
angles = [0, theta_d, pi/2, pi-theta_d, pi, -pi+theta_d, -pi/2, -theta_d];
direction_indices = struct('east', 1, 'northeast', 2, 'north', 3, ...
                           'northwest', 4, 'west', 5, 'southwest', 6, ...
                           'south', 7, 'southeast', 8);


inward_direction_index = mod(direction_indices.(direction) + 3, 8) + 1;
inward_direction_index_plus_1 = mod(inward_direction_index, 8) + 1;
inward_direction_index_minus_1 = mod(inward_direction_index - 2, 8) + 1;

inward_angle = angles(inward_direction_index);
next_inward_angle = angles(inward_direction_index_plus_1);
prev_inward_angle = angles(inward_direction_index_minus_1);

interp_table_x = [-angular_difference(inward_angle, prev_inward_angle), ...
                  0, ...
                  angular_difference(next_inward_angle, inward_angle)];
interp_table_y = [0 1 0];

w = interp1(interp_table_x, interp_table_y, angular_difference(R, inward_angle));
w(isnan(w)) = 0;
%--------------------------------------------------------------------------

%==========================================================================
function d = angular_difference(theta1, theta2)
%
% Given arrays of angles in radians, calculate the absolute angular differences
% between the values in theta1 and the values in theta2.  Return the difference
% values in the half-open interval [-pi, pi).
%
% Remove multiples of 2*pi in the difference calculation.  In other words,
% angular_difference(0,4*pi) returns 0, and angular_difference(0,2*pi + pi/4)
% returns -pi/4.

d = mod(theta1 - theta2 + pi, 2*pi) - pi;
%--------------------------------------------------------------------------
