function [mdl,int,intci,predstats,rank] = bayesloglinear(P,c,varargin)

%BAYESLOGLINEAR Bayesian analysis of a loglinear point process model 
%
% Syntax
%
%     [mdl,int,intci,predstats] = bayesloglinear(P,c,pn,pv,...)
%
% Description
%
%     Loglinear models embrace numerous models that can be fitted to
%     homogeneous and inhomogeneous Poisson processes on river networks. 
%     bayesloglinear uses bayesreg, a toolbox for fitting and evaluating 
%     Bayesian penalized regression models with continuous shrinkage prior 
%     densities.
%
% Input arguments
%
%     P      instance of PPS
%     c      covariates (node-attribute lists)
%     
%     Parameter name/value pairs
%     
%     'verbose'     true
%     'prior'       {'horseshoe'}. Other options are 'g','ridge','lasso',
%                   'horseshoe+' 
%     'predci'      [5 95]. prediction confidence intervals
%     'nsamples'    1000. number of posterior samples
%     'thin'        5. level of thinning
%     'varnames'    {}. cell array of variable names
%     'sortrank'    false. display variables in rank order
%
% Output arguments
%
%     mdl    model (GeneralizedLinearModel)
%     int    node-attribute list with modelled intensities
%     intc   credible intervals (5 and 95%) of intensities
%     predstats see bayesreg
%     
% Example   
% 
% See also: PPS, PPS/random, PPS/fitloglinear, bayesreg
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. February, 2019

% Check, whether bayesreg is available
if exist('bayesreg','file') ~= 2
    url = 'https://in.mathworks.com/matlabcentral/fileexchange/60823-bayesian-penalized-regression-with-continuous-shrinkage-prio';
    error('TopoToolbox:bayesloglinear',...
        ['BayesReg Toolbox must be on the search path. The toolbox\n'...
         'can be downloaded from the <a href="' url '">File Exchange</a>.'])
end

p = inputParser;
p.KeepUnmatched = true;
addParameter(p,'distribution','binomial');
addParameter(p,'verbose',true);
addParameter(p,'prior','horseshoe');
addParameter(p,'predci',[5 95],@(x) (numel(x) == 2));
addParameter(p,'nsamples',1000)
addParameter(p,'thin', 5)
addParameter(p,'varnames',{})
addParameter(p,'sortrank',false)
parse(p,varargin{:});

% Options for fitglm and stepwise
bayesregopts = p.Unmatched;
bayesregopts = expandstruct(bayesregopts);

% Get covariates
X = getcovariate(P,c);
% Get response variable
y = zeros(numel(P.S.x),1);
y(P.PP) = 1;

% [b,FitInfo] = lassoglm(zscore(X),y,'binomial');

[mdl.beta,mdl.beta0,mdl.retval] = ...
    bayesreg(X,y,p.Results.distribution,p.Results.prior,...
                    'display',p.Results.verbose,...
                    'nsamples',p.Results.nsamples,...
                    'thin',p.Results.thin,...
                    'varnames',p.Results.varnames,...
                    'sortrank',p.Results.sortrank,...
                    bayesregopts{:});
                
[pred, predstats] = br_predict(X, mdl.beta, mdl.beta0, mdl.retval, ...
                    'CI',p.Results.predci,...
                    'ytest', y, ...
                    'display',p.Results.verbose);
                
rank = bfr(mdl.beta);

% modelled intensities
d   = distance(P.S,'node_to_node');
d   = mean(d);
int = pred.prob_1./d;
intci = [pred.prob_1_CI5_./d pred.prob_1_CI95_./d];


end

function pnpv = expandstruct(s)

pn = fieldnames(s);
pv = struct2cell(s);

pnpv = [pn pv];
pnpv = pnpv';
pnpv = pnpv(:)';
end