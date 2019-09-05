function [mn,results] = mnoptim(S,DEM,A,varargin)

%MNOPTIM Bayesian optimization of the mn ratio
%
% Syntax
%   
%     [mn,results] = mnoptim(S,DEM,A)
%     [mn,results] = mnoptim(S,z,a)
%     [mn,results] = mnoptim(...,pn,pv,...)
%
% Description
%
%     mnoptim uses Bayesian optimization to find an optimal mn-ratio by
%     minimizing a cross-validation loss. Cross-validation is performed
%     using repeated 2-fold validation on individual stream networks
%
% Input arguments
%
%     S      STREAMobj
%     DEM    Digital elevation model (GRIDobj)
%     A      upstream area (in pixels) as returned by the function flowacc
%            (GRIDobj)
%     z      node-attribute list of elevations
%     a      node-attribute list of upstream areas (nr of upstream pixels)
%
% Parameter-name value pairs
%
%     'optvar'    {'mn'},'m','n','m&n'
%                 The variable to be optimized. 'mn' optimizes the mn ratio,
%                 'm' finds an optimal value of m for a given value of n.
%                 'n' finds an optimal value of n for a given value of m.
%                 'm&n' find an optimal value of m and n.
%     'n'         value of n, only applicable if 'optvar' is 'm'. Default
%                 value is 1.
%     'm'         value of m, only applicable if 'optvar' is 'n'. Default
%                 value is 0.5. 
%     'nrange'    search range of 'n'. Only applicable if 'optvar' is 'n' or
%                 'm&n'. Default is [0.8 1.5].
%     'mrange'    search range of 'm'. Only applicable if 'optvar' is 'm'.
%                 'm&n'.Default is [0.3 1].
%     'mnrange'   search range of 'mn'. Only applicable if 'optvar' is 'm'.
%                 Default is [.1 1].
%     'a0'        reference area (see chitransform). Default is 1e6 m^2.
%     'lossfun'   function to be minimized. Default is 'rmse'. Others are
%                 'linear' and 'coeffdeterm'.
%     'UseParallel' {true} or false.
%     'crossval'  {true} or false.   
%
% Other bayesopt parameter name/value pairs (see help bayesopt for details)
%
%     'MaxObjectiveEvaluations'
%
%
% Output arguments
%
%     mn       table with results (best point)
%     results  BayesianOptimization object. Continue with optimization
%              with the command resultsnew = resume(results).
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM);
%     A  = flowacc(FD);
%     S  = STREAMobj(FD,'minarea',1e6,'unit','map');
%     S  = removeedgeeffects(S,FD,DEM);
%     [mn,results] = mnoptim(S,DEM,A,'optvar','mn','crossval',true);
%     
% See also: chiplot, chitransform, slopearea
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 24. October, 2017

p = inputParser;         
p.FunctionName = 'STREAMobj/mnoptim';
addRequired(p,'S',@(x) isa(x,'STREAMobj'));
addRequired(p,'DEM',@(x)  isa(x,'GRIDobj') || isnal(S,x));
addRequired(p,'A', @(x) isa(x,'GRIDobj') || isnal(S,x));

addParamValue(p,'mnrange',[.1 1],@(x) numel(x)==2);
addParamValue(p,'optvar','mn');
addParamValue(p,'n',1);
addParamValue(p,'m',0.5);
addParamValue(p,'nrange',[0.8 1.5]);
addParamValue(p,'mrange',[0.3 1]);
addParamValue(p,'a0',1e6,@(x) isscalar(x) && isnumeric(x));
addParamValue(p,'lossfun','rmse',@(x) ischar(x) || isa(x,'function_handle'));
addParamValue(p,'plot',false);
addParamValue(p,'crossval',true);
addParamValue(p,'verbose',0);
addParamValue(p,'UseParallel',true);
addParamValue(p,'MaxObjectiveEvaluations',30);



parse(p,S,DEM,A,varargin{:});
a0 = p.Results.a0;

% get node attribute list with elevation values
if isa(DEM,'GRIDobj')
    validatealignment(S,DEM);
    z = getnal(S,DEM);
elseif isnal(S,DEM)
    z = DEM;
else
    error('Imcompatible format of second input argument')
end

% get node attribute list with flow accumulation values
if isa(A,'GRIDobj')
    validatealignment(S,A);
    a = getnal(S,A);
elseif isnal(S,A)
    a = A;
else
    error('Imcompatible format of second input argument')
end
clear A DEM

% set the base level of all streams to zero
zbase = z;
for r = numel(S.ixc):-1:1
    zbase(S.ix(r)) = zbase(S.ixc(r));
end
z = z-zbase;

%% ----- beta ------------------------------------------------------------
% set all rivers to have the same gradient 





%% -----------------------------------------------------------------------

if p.Results.crossval
    [~,locb] = STREAMobj2cell(S);
    cv = true;
    % Number of connected components
    nrcc = numel(locb);
    if nrcc == 1
        cv = false;
        error('TopoToolbox:mnoptim',['Cross-validation not possible. There is only \n' ...
                                       'one connected component in the stream network.\n' ...
                                       'Set option ''crossval'' to false.']);
    end
else
    cv = false;
end

% get lossfunction
if ischar(p.Results.lossfun)
    switch p.Results.lossfun
        case 'rmse'
            lossfun = @(x,xhat)mean(sum((x-xhat).^2));
        case 'linear'
            lossfun = @(x,xhat)var((x+1)./(xhat+1));
        case 'linearcc'
            cc = conncomps(S);
            lossfun = @(x,xhat) sum(accumarray(cc,(x+1)./(xhat+1),[],@std)).^2;
%             lossfun = @(x,xhat)var((x+1)./(xhat+1));
            
        case 'coeffdeterm'
            take2elem = @(x) x(2);
            lossfun = @(x,xhat)1-(take2elem(corrcoef(x,xhat))).^2;
        otherwise
            error('unknown loss function')
    end
else
    lossfun = p.Results.lossfun;
end

% Bayes
optvar = p.Results.optvar;
isdeterm = ~cv;
pp   = gcp('nocreate');
opts = {'IsObjectiveDeterministic', isdeterm, ...
        'Verbose', p.Results.verbose,...
        'UseParallel',~(isempty(pp) & ~p.Results.UseParallel),...
        'MaxObjectiveEvaluations',p.Results.MaxObjectiveEvaluations,...
        };
        
switch optvar
    case 'mn'
        mn      = optimizableVariable('mn',p.Results.mnrange);    
        results = bayesopt(@(x) fun(x),mn,opts{:});
    case 'm'
        m       = optimizableVariable('m',p.Results.mrange);
        n       = p.Results.n;
        results = bayesopt(@(x) fun(x),m,opts{:});    
    case 'n'
        n       = optimizableVariable('n',p.Results.nrange);
        m       = p.Results.m;
        results = bayesopt(@(x) fun(x),n,opts{:});
    case 'm&n'
        m       = optimizableVariable('m',p.Results.mrange);
        n       = optimizableVariable('n',p.Results.nrange);
        results = bayesopt(@(x) fun(x),[m n],opts{:});
               
end

mn      = bestPoint(results);



% loss function
function lss = fun(x)


switch optvar
    case 'mn'
        c = chitransform(S,a,'mn',x.mn,'a0',a0);    
    case 'm'
        c = chitransform(S,a,'mn',x.m/n,'a0',a0);  
    case 'n'
        c = chitransform(S,a,'mn',m/x.n,'a0',a0);  
    case 'm&n'
        c = chitransform(S,a,'mn',x.m/x.n,'a0',a0);  
end

if ~cv    
    % no cross-validation
    b = z\c;
    zhat = b*c;
    lss = lossfun(z/max(z),zhat/max(zhat));
else
    % cross validation
    LSS = zeros(5,1);
    for iter = 1:numel(LSS)
        conncomptrain = randperm(nrcc,ceil(nrcc/2));
        CC            = false(nrcc,1);
        CC(conncomptrain) = true;
        ixtrain       = vertcat(locb{CC});
        b = z(ixtrain)\c(ixtrain);
        ixval   = vertcat(locb{~CC});
        zhat = b*c(ixval);
        zt   = z(ixval);
        LSS(iter)  = lossfun(zt/max(zt),zhat/max(zhat));
%         LSS(iter)  = lossfun(zt/prctile(zt,75),zhat/prctile(zhat,75));
    end
    lss = mean(LSS);
end
    

end
end

