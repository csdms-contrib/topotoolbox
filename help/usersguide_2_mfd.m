%% TopoToolbox 2 - calculating the multiple flow direction
%
% <<TopoToolbox.png>>
%
%
% With the introduction of TopoToolbox 2, the new FLOWobj does only support
% the single (D8) flow direction algorithm. If you have been a user of
% TopoToolbox 1 and require the multiple flow direction (MFD) algorithm,
% please read these pages.

%% Requirements
%
% If you want to use the MFD algorithm with TopoToolbox 2, you must have a 
% previous version (e.g. TopoToolbox 1.6) on your search path.

%% Here is how to do it
%

% get the coordinate vectors and the DEM in the format required by the
% function flowdir
[dem,X,Y] = GRIDobj2mat(fillsinks(DEM));
% use meshgrid to calculate coordinate matrices
[X,Y] = meshgrid(X,Y);
% calculate the flow direction matrix M
M = flowdir(X,Y,double(dem));

%%
% Now you should calculate the desired topographic attribute directly from
% M, e.g. flow accumulation

A  = flowacc(M,size(dem));

%%
% and convert the output back to a GRIDobj.
A  = GRIDobj(X,Y,A);
% if there is projection information in the property .georef, you can copy
% it from the DEM.
A.georef = DEM.georef;
% and plot it
imageschs(DEM,log(A));
%%
% and this is it.

%% Future developments
% The fact, that TT2 currently does not support the MFD algorithm is
% unsatisfactory. We are currently working on including MFD in
% future version of TT.



