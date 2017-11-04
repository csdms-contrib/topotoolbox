function DEM = quantcarve(FD,DEM,tau)

%QUANTCARVE quantile carving
%
% Syntax
%
%     DEMc = quantcarve(S,DEM,tau)
%
% Description
%
%     Elevation values along stream networks are frequently affected by
%     large scatter, often as a result of data artifacts or errors. This
%     function returns a node attribute list of elevations calculated by
%     carving the DEM. Conversely to conventional carving, quantcarve will
%     not run along minimas of the DEM. Instead, quantcarve returns a
%     profile that runs along the tau's quantile of elevation conditional 
%     horizontal distance of the river profile.
%
%     The function uses linprog from the Optimization Toolbox.
%
%     Applying quantile carving to entire DEMs is computationally
%     intense. Although this function will be run in parallel, it might
%     take a long while to execute, sometimes it will even throw an error.
%     In this case, consider to derive a stream network (STREAMobj) and
%     apply quantile carving to the stream network only (see
%     STREAMobj/quantcarve). Then map the carved values back to the DEM.
%
% Input parameters
%
%     FD     FLOWobj
%     DEM    Digital elevation model (GRIDobj)
%     tau    quantile (default is 0.5)
%
% Output parameters
%
%     DEMc   carved DEM
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     DEMc = quantcarve(FD,DEM,0.1);
%     imageschs(DEMc,DEM-DEMc)
%
%
% See also: STREAMobj/quantcarve, FLOWobj/imposemin, STREAMobj/crs
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 18. July, 2017


% check input arguments
narginchk(2,inf)
if nargin == 2
    tau = 0.5;
end
validateattributes(tau,{'numeric'},{'>',0,'<',1},'FLOWobj/quantcarve','tau',3);

% identify filled pixels
I = fillsinks(DEM) > DEM;
I = dilate(I,true(3)) & ~isnan(DEM);

% derive STREAMobj with channelheads upstream of filled pixels
I = influencemap(FD,I);
S = STREAMobj(FD,I);

% provide some output
disp(['The total number of nodes is ' num2str(numel(S.IXgrid))])

% quantile carving
zs = quantcarve(S,DEM,tau,'split',2);

% map values back to DEM
DEM.Z(S.IXgrid) = cast(zs,class(DEM.Z));