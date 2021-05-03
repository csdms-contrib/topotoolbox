function [a,mask] = maplateral(S,A,dist,aggfun,varargin)

%MAPLATERAL map values of regions adjacent to streams to stream network
%
% Syntax
%
%     a = maplateral(S,A,dist,aggfun)
%     [a,mask] = ...     
%     a = maplateral(...,pn,pv,...)
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
%     dist    scalar indicating maximum distance from stream network or
%             two element vector with minimum and maximum distance or
%             node-attribute list with maximum distance values for each
%             node in the river network or
%             two columns of node-attribute lists indicating
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
%     flat              {false} or true. If true, than the ends of the
%                       buffer at channel heads or outlets are flat,
%                       otherwise they are round.
%     fillval           fill value if stream pixel has no nearest neighbors
%                       By default, this is nan.
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
% Date: 23. February, 2021

% Input checking and parsing
validatealignment(S,A)
p = inputParser;
p.FunctionName = 'STREAMobj/maplateral';
addParamValue(p,'excludestream',true,@(x) isscalar(x));
addParamValue(p,'inpaintnans',true,@(x) isscalar(x));
addParamValue(p,'flat',false,@(x) isscalar(x));
addParamValue(p,'fillval',nan);
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

% Create mask
if isscalar(dist)
    I  = D<=dist;
elseif numel(dist) == 2
    dist = sort(dist);
    I  = D<=dist(2) & D>=dist(1);
elseif isnal(S,dist(:,1))
    if size(dist,2) == 1
        W = STREAMobj2GRIDobj(S,dist);
        I = D <= W.Z(L);
    else
        dist = sort(dist,2,'ascend');
        W1 = STREAMobj2GRIDobj(S,dist(:,1));
        W2 = STREAMobj2GRIDobj(S,dist(:,2));
        I  = D <= W2.Z(L) & D >= W1.Z(L);
    end
else
    error('dist must be a scalar, a two-element vector, or a node-attribute list');
end

% flat tops?
if p.Results.flat
    endpoints = streampoi(S,{'channelhead','outlet'},'logical');
    II = ismember(L,S.IXgrid(endpoints));
    I(II) = false;
    I  = (imdilate(I,ones(3)) & II) | I;
end

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
        a(:,r)  = accumarray(locb,double(A),size(S.IXgrid),aggfun{r},p.Results.fillval);
    end
else
    a  = accumarray(locb,double(A),size(S.IXgrid),aggfun,p.Results.fillval);
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

