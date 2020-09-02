function [mn,cm,zm,zsd] = mnoptimvar(S,DEM,A,varargin)

%MNOPTIMVAR optimize the mn ratio using minimum variance method
%
% Syntax
%
%     mn = mnoptimvar(S,DEM,A)
%     mn = mnoptimvar(S,z,a)
%     mn = mnoptimvar(...,pn,pv,...)
%     [mn,cm,zm,zsd] = ...
%
% Description
%
%     mnoptimvar finds an optimal value of the mn-ratio by minimizing the
%     variance of elevation values conditional on chi. The procedure places
%     nodes of the river network into chi-distance bins and calculates the
%     variance of elevation values in each bin. The objective function is
%     the weighted mean of these variances.
%
% Input arguments
%
%     S     STREAMobj
%     DEM   Digital elevation model
%     A     Flow accumulation grid 
%     z     node-attribute list of elevations
%     a     node-attribute list of flow accumulation
%
%     Parameter name/value pairs
%
%     'varfun'   function to calculate elevation variability in chi-distance 
%                bins. Default is the interquartile range (@iqr), but any 
%                other measure of variability may be suitable, e.g. @var,
%                @std, @mad, @robustcov, @range.
%     'mn0'      starting value for mn {0.5}.
%     'distbins' nr of distance bins {100}.
%     'minstream' minimum number of streams that should be included in
%                calculation a variance value {2}.
%     'funoptim' {'fminsearch'} (so far only choice)
%     'optimize' {true} or false. If false, then the function won't perform
%                an optimization, but will use mn0 as mn-ratio.
%     'plot'     {true} or false
%     'zerobaselevel' {false} or true. Set all outlet elevations to zero.
%     'a0'       reference area {1e6}
%     
% Output arguments
%
%     mn    mn-ratio
%     cm    binned chi values
%     zm    average elevation in chi bins
%     zsd   standard deviation of elevations in chi bins
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     A        = flowacc(FD);
%     S   = STREAMobj(FD,'minarea',1000);
%     S   = klargestconncomps(S);
%     mn  = mnoptimvar(S,DEM,A);
%     
%
% See also: STREAMobj, STREAMobj/mnoptim
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 29. January, 2020

% parse inputs
p = inputParser;
addParameter(p,'varfun',@iqr,@(x) isa(x,'function_handle'));
addParameter(p,'distbins',100,@(x) isscalar(x) && x > 1);
addParameter(p,'minstreams',2,@(x) isscalar(x) && x > 1);
addParameter(p,'mn0',0.5,@(x) isscalar(x) && x > 0);
addParameter(p,'funoptim','fminsearch')
addParameter(p,'plot',true)
addParameter(p,'optimize',true);
addParameter(p,'zerobaselevel',false);
addParameter(p,'a0',1e6);
parse(p,varargin{:});

% get parameters
varfun = p.Results.varfun;
distbins = p.Results.distbins;
minstreams = p.Results.minstreams;
a0 = p.Results.a0;
    

% get node attribute list with flow accumulation values
if isa(A,'GRIDobj')
    validatealignment(S,A);
    a = getnal(S,A);
elseif isnal(S,A)
    a = A;
else
    error('Imcompatible format of second input argument')
end

% get node attribute list with elevation values
if isa(DEM,'GRIDobj')
    validatealignment(S,DEM);
    z = getnal(S,DEM);
elseif isnal(S,DEM)
    z = DEM;
else
    error('Imcompatible format of second input argument')
end

% z must be double precision
z = double(z);
if p.Results.zerobaselevel
    z = zerobaselevel(S,z);
end
copyz = z;
z = z - min(z) + 1;

if p.Results.optimize
% calculate stream 
label = labelreach(S);

% choose optimizer
switch lower(p.Results.funoptim)
    case 'fminsearch'
        mn = fminsearch(@(mn) chivar(mn), log(p.Results.mn0));
        mn = exp(mn);
end
else 
    mn = p.Results.mn0;
end

if nargout > 1 || p.Results.plot
    c   = chitransform(S,a,'mn',mn,'a0',a0);
    [~,~,bin] = histcounts(c,distbins);
    zm  = accumarray(bin,copyz,[],@mean,0);
    zsd = accumarray(bin,copyz,[],@std,0);
    cm  = accumarray(bin,c,[],@mean,0);
    
end

if p.Results.plot
    plotdz(S,copyz,'distance',c,'color',[.7 .7 .7])
    hold on
    errorbar(cm,zm,zsd,'k.')
    plot(cm,zm,'ks-');
    hold off
    xlabel('\chi [m]')
end
    

%% objective function
function v = chivar(mn)
    
    % calculate chitransform
    c = chitransform(S,a,'mn',exp(mn),'a0',a0);
    % bin c values into distance bins.
    [~,~,bin] = histcounts(c,distbins);
    
    % calculate the average values within each bin and tributary label
    melev = accumarray([label bin],z,[],@mean,0,true);
    
    % calculate cross-tributary variability using varfun
    n     = sum(spones(melev))';
    [~,j,melev] = find(melev);
    v    = accumarray(j(n(j)>=minstreams),melev(n(j)>=minstreams),...
                      [numel(n) 1],varfun,nan);
    
    % exclude those variances that have  less than minstreams values              
    I    = n < minstreams;
    v(I) = [];
    n(I) = [];
    
    % calculate the weighted average of all variances. The weighting is
    % according to the number of streams used to calculate the variance in
    % each bin.
    n = n/sum(n);
    v    = n'*v;
end
end
