function zs = aggregate(S,DEM,varargin)

%AGGREGATE aggregating values of reaches
%
% Syntax
%
%     as = aggregate(S,A)
%     as = aggregate(S,a)
%     as = aggregate(...,pn,pv,...)
%
% Description
%
%     Values along stream networks are frequently affected by scatter.
%     This function removes scatter by averaging over values in reaches
%     between confluences or reaches of equal length.
%
% Input parameters
%
%     S        STREAMobj
%     A        GRIDobj
%     a        node attribute list
%
% Parameter name/value pairs
%
%     'method'     {'reach'}, 'betweenconfluences','drainagebasins'
%     'split'      {false} or true. True will identify individual drainage
%                  basins and process each individually in parallel (requires
%                  the parallel processing toolbox).
%     'seglength'  segment length (default 10*S.cellsize)
%     'aggfun'     anonymous function as aggregation function. Default is
%                  @mean. The function must take a vector and return a scalar
%                  (e.g. @max, @min, @std, @(x) prctile(x,25), ...)
%
% Output parameters
%
%     as     node attribute list with aggregated values
%
% Example: Calculating ksn values
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S);
%     g = gradient(S,imposemin(S,DEM));
%     A = flowacc(FD);
%     a = getnal(S,A)*(DEM.cellsize^2);
%     ksn = g./(a.^(-0.45));
%     ksn = aggregate(S,ksn,'seglength',2000);
%     plotc(S,ksn)
%     axis image
%     colormap(jet)
%     colorbar
%
% See also: STREAMobj/labelreach, STREAMobj/smooth
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 19. June, 2017


% check and parse inputs
narginchk(2,inf)

p = inputParser;
p.FunctionName = 'STREAMobj/aggregate';
addParameter(p,'method','reach'); 
addParameter(p,'split',false);
% parameters for block and reach
addParameter(p,'seglength',S.cellsize*11);
addParameter(p,'aggfun',@mean);

parse(p,varargin{:});
method = validatestring(p.Results.method,{'betweenconfluences','reach','drainagebasins'});

% get node attribute list with elevation values
if isa(DEM,'GRIDobj')
    validatealignment(S,DEM);
    z = getnal(S,DEM);
elseif isnal(S,DEM)
    z = DEM;
else
    error('Imcompatible format of second input argument')
end

% run in parallel if wanted
if p.Results.split
    params       = p.Results;
    params.split = false;
    [CS,locS]    = STREAMobj2cell(S);
    Cz           = cellfun(@(ix) z(ix),locS,'UniformOutput',false);
    Czs          = cell(size(CS));
    parfor r = 1:numel(CS)
        Czs{r} = aggregate(CS{r},Cz{r},params);
    end
    
    zs = nan(size(z));
    for r = 1:numel(CS)
        zs(locS{r}) = Czs{r};
    end
    return
end

%% Aggregating starts here

switch method
    case 'betweenconfluences'
        S2 = split(S);
        [c,n] = conncomps(S2);
        za = accumarray(c,z,[n 1],p.Results.aggfun,nan,false);
        zs = za(c);
    case 'reach'
        label = labelreach(S,'seglength',p.Results.seglength);   
        zm    = accumarray(label,z,[max(label) 1],p.Results.aggfun,nan);
        zs    = zm(label);
    case 'drainagebasins'
        label = conncomps(S);
        zm    = accumarray(label,z,[max(label) 1],p.Results.aggfun,nan);
        zs    = zm(label);
end

end
        
