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
% Algorithm
%
%     This algorithm uses quantile carving to smooth the data. The
%     algorithm is described in Schwanghart and Scherler (2017) (Eq. A12).
%     
% References
%
%     Schwanghart, W., Scherler, D., 2017. Bumps in river profiles: 
%     uncertainty assessment and smoothing using quantile regression 
%     techniques. Earth Surface Dynamics, 5, 821-839. 
%     [DOI: 10.5194/esurf-5-821-2017]
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
I = influencemap(FD,I);

% I = GRIDobj(DEM,'logical');
% I.Z(:,:) = true;

% derive STREAMobj with channelheads upstream of filled pixels
CFD = FLOWobj2cell(FD);
hassinks   = cellfun(@(x) any(I.Z(x.ix)),CFD);
CFD = CFD(hassinks);

nd = numel(CFD);
CFD = CFD(randperm(nd));
ndstring = num2str(nd);
h  = waitbar(0,['Please wait... (0/' ndstring ')']);

for r = 1:nd
    rstring = num2str(r);
    waitbar(r/nd,h,['Calculate stream network (' rstring '/' ndstring ')']);
    S = STREAMobj(CFD{r},I);

    % quantile carving
    waitbar(r/nd,h,['Start quantile carving (nrnodes = ' num2str(numel(S.IXgrid))  ', ' rstring '/' ndstring ')']);
    zs = quantcarve(S,DEM,tau,'split',0,'waitbar',false);
    waitbar(r/nd,h,['Finished quantile carving (' rstring '/' ndstring ')']);
    % map values back to DEM
    DEM.Z(S.IXgrid) = cast(zs,class(DEM.Z));
end
close(h);