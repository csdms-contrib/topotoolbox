function mask = border_nans(E)
%border_nans Find NaNs connected to the DEM border
%
%   mask = border_nans(E) finds all connected NaN components that touch the
%   border of E, the digital elevation model (DEM). mask is a logical
%   matrix the same size as E.  True values in mask correspond to border
%   NaN locations. 
%
%   Example
%   -------
%      E = magic(5);
%      E(2:5,1:2) = NaN;
%      E(3:4,4) = NaN
%      border_nans(E)
%
%   See also pixel_flow, upslope_area.

%   Steven L. Eddins
%   $Revision: 1.1 $  $Date: 2007/10/02 15:49:17 $
%   Copyright 2007 The MathWorks, Inc.

nan_mask = isnan(E);
mask = nan_mask & ~imclearborder(nan_mask);
