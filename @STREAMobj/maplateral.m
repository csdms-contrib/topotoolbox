function [a,mask] = maplateral(S,A,dist,aggfun,varargin)

%MAPLATERAL map values of regions adjacent to streams to stream network
%
% Syntax
%
%     a = maplateral(S,A,dist,aggfun)
%     [a,mask] = ...     
%
% Description
%
%     Similarly to swath profiles, maplateral maps nearest values from a  
%     grid A to a stream network S within a distance dist and using the 
%     aggregation function aggfun. The function returns a node attribute 
%     list a.
%
% Input arguments
%
%     S       STREAMobj
%     A       GRIDobj
%     dist    scalar indicating maximum distance from stream network.
%     aggfun  anonymous aggregation function (e.g. @mean). aggfun can be a
%             cell array of anonymous functions (e.g. {@min @max}). In this
%             case, the output a will have as many columns as there are
%             functions in aggfun.
%
%     Parameter name/value pairs
%
%     excludestream     {true} or false. If true, grid values of the 
%                       stream network are excluded.
%     inpaintnans       {true} or false. If true, missing values in the
%                       resulting nal will be interpolated (see STREAMobj/
%                       inpaintnans).
%
% Output arguments
%
%     a       node-attribute list
%     mask    logical GRIDobj with true values referring to calculating 
%             aggregated values. 
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S);
%     S = trunk(S);
%     g = maplateral(S,gradient8(DEM),100,{@min @max @mean});
%     plotdz(S,DEM)
%     yyaxis right
%     plotdzshaded(S,g(:,[1 2])
%     hold on
%     plotdz(S,g(:,3))
%     ylabel('Hillslope gradient [-]')
%
% See also: STREAMobj, SWATHobj, STREAMobj/smooth
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 12. June, 2018

% Input checking and parsing
validatealignment(S,A)
p = inputParser;
p.FunctionName = 'STREAMobj/maplateral';
addParamValue(p,'excludestream',true,@(x) isscalar(x));
addParamValue(p,'inpaintnans',true,@(x) isscalar(x));
parse(p,varargin{:});

% create mask, if required
if nargout == 2
    mask = A;
end

% Stream raster
B     = false(A.size);
B(S.IXgrid) = true;
% Eucl. distance transform
[D,L] = bwdist(B,'e');
dist  = ceil(dist/S.cellsize);

if isscalar(dist)
    dist = [0 dist];
end

% mask
I  = D<=dist(2) & D>=dist(1);

if p.Results.excludestream
    I(S.IXgrid) = false;
end

% grid values
A  = A.Z(I);
% closest pixel map
L  = L(I);

A  = A(:);
L  = L(:);

% aggregate values at each stream node
[~,locb] = ismember(L,S.IXgrid);
if iscell(aggfun)
    a = nan(numel(S.IXgrid),numel(aggfun));
    for r = 1:numel(aggfun)
        a(:,r)  = accumarray(locb,double(A),size(S.IXgrid),aggfun{r},nan);
    end
else
    a  = accumarray(locb,double(A),size(S.IXgrid),aggfun,nan);
end

% inpaint nans
if p.Results.inpaintnans
    for r = 1:size(a,2)
        a(:,r) = inpaintnans(S,a(:,r));
    end
end

% second output if required
if nargout == 2
    mask.Z = I;
    mask.name = 'mask';
end

