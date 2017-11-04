Version 1.4.2
=============
16-Sep-2009

* Updated unit tests to use assertElementsAlmostEqual.

Version 1.4.1
=============
10-Mar-2009

* Fixed M-file help typo in dependence_map.m.

* Allowed input matrix E to contain NaNs in postprocess_plateaus.m.

* Added unit test directory.

Version 1.4
===========
25-Sep-2008

* Modified fill_sinks.m to call imfill.m with 8-connectivity
  instead of 4-connectivity.

Version 1.3
===========
14-Feb-2008

* Improved handling of groups of NaNs that touch the DEM border
  so that the dependence_map and influence_map calculations are
  correct.  Now flow_matrix is computed so that border NaN pixels
  have zero flow weights to and from all their neighbors.  As a
  nice side effect of the change, flow_matrix.m is now faster
  for datasets that have border NaN pixels.

Version 1.2
===========
02-Oct-2007

* Changed handling of groups of NaNs that touch the DEM border.

* Added border_nans function.

Version 1.1
===========
06-Aug-2007

* Incompatible change made to upslope_area.  This function no longer
  "flattens" the upslope areas computed for plateaus.

* New function: postprocess_plateaus.  This function flattens the upslope
  areas computed for plateaus.  This function was formerly a part of
  upslope_area.m


Version 1.0
===========
02-Aug-2007
Initial release
