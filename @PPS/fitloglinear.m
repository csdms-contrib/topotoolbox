function [mdl,int] = fitloglinear(P,c,varargin)

%FITLOGLINEAR fit loglinear model to point pattern
%
% Syntax
%
%     [mdl,int] = fitloglinear(P,c)
%     [mdl,int] = fitloglinear(P,c,pn,pv,...)
%
% Description
%
%     Loglinear models embrace numerous models that can be fitted to
%     homogeneous and inhomogeneous Poisson processes on river networks. 
%
% Input arguments
%
%     P      instance of PPS
%     c      covariates (node-attribute lists)
%     
%     Parameter name/value pairs
%     
%     'stepwise'   {false} or true. If true, fitloglinear uses stepwiseglm
%                  to fit the model.
%     'modelspec'  see fitglm
%     
%     In addition, fitloglinear accepts parameter name/value pairs of
%     fitglm or stepwiseglm.
%
% Output arguments
%
%     mdl    model (GeneralizedLinearModel)
%     int    node-attribute list with modelled intensities
%     
% Example   
% 
% See also: PPS, PPS/random
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. February, 2020

% See also: fitglm, stepwiseglm

p = inputParser;
p.KeepUnmatched = true;
addParameter(p,'stepwise',false,@(x) isscalar(x));
addParameter(p,'modelspec','linear');
addParameter(p,'distribution','binomial');
parse(p,varargin{:});

% Get covariates
X = getcovariate(P,c);
% Get response variable
y = zeros(numel(P.S.x),1);
y(P.PP) = 1;

% Options for fitglm and stepwise
glmopts = p.Unmatched;
glmopts = expandstruct(glmopts);

% Covariates can be supplied as table. This needs some handling.
if istable(X)
    tbl = [X table(y,'variablenames',{'response'})];
    inp = {tbl};
else
    inp = {X y};
end


if ~p.Results.stepwise
    
    mdl = fitglm(inp{:},p.Results.modelspec,...
        'Distribution',p.Results.distribution,glmopts{:});
    
else

    % Use stepwise GLM
    mdl = stepwiseglm(inp{:},p.Results.modelspec,...
        'Distribution',p.Results.distribution,glmopts{:});
    
end

% modelled intensities
p   = predict(mdl,X);
d   = distance(P.S,'node_to_node');
d   = mean(d);
int = p./d; %.cellsize;

end

function pnpv = expandstruct(s)

pn = fieldnames(s);
pv = struct2cell(s);

pnpv = [pn pv];
pnpv = pnpv';
pnpv = pnpv(:)';
end