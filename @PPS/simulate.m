function CP = simulate(P,varargin)

%SIMULATE simulate point pattern using random thinning
%
% Syntax
%
%     PSIM = simulate(P,'intensity',int)
%     PSIM = simulate(P,'prob',p)
%     C    = simulate(P,...,'nsim',n)
%
% Description
%
%     
%
% Input arguments
%
%     P      point pattern (PPS)
%
%     Parameter name/value pairs
%
%     'intensity'    node-attribute list with intensities or scalar
%                    intensity
%     'prob'         probability of a point to remain. Thins an existing
%                    point pattern based on a scalar probability or a
%                    probability for each point in P.
%     'nsim'         number of simulations
%
% Algorithm
%
%     simulate with the option 'intensity' randomly generates a point
%     pattern with Poisson distribution and then uses random thinning to
%     remove points. Note that very high intensities can cause problems.
%
% See also: PPS, PPS/convhull
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. February, 2019

narginchk(3,inf)

% Check input arguments
p = inputParser;
p.FunctionName = 'PPS/simulate';
addRequired(p,'P',@(x) isa(x,'PPS'));
addParameter(p,'intensity',[],@(x) isnal(P.S,x) || isscalar(x) || isempty(x))
addParameter(p,'prob',[],@(x) (numel(x) == npoints(P) || isscalar(x)) && all(x>=0 & x<=1) || isempty(x));
addParameter(p,'nsim',1);

% Parse
parse(p,P,varargin{:});

% run several simulations
if p.Results.nsim > 1 
    CP   = cell(1,p.Results.nsim);
    for r = 1:p.Results.nsim
        CP{r} = simulate(P,'intensity',p.Results.intensity,'nsim',1,'prob',p.Results.prob);
    end
    return
end

% intensity
int = p.Results.intensity;
if ~isempty(int)
    if isscalar(int)
        int = repmat(int,numel(P.S.IXgrid),1);
    end
end

if ~isempty(int)
    
    % simulate a poisson process with the maximum intensity in int
    maxint = max(int);
    PN     = PPS(P.S,'rpois',maxint);
    P.PP   = PN.PP;
    
    % random thinning
    % The probability for a point to remain in the list is the
    % ratio between intensity and the maximum intensity
    p_remain = int(P.PP)./maxint;
    
    
    
else
    
    if isscalar(p.Results.prob)
        p_remain = repmat(p.Results.prob,npoints(P),1);
    else
        p_remain = p.Results.prob;
    end
    
end

remain   = binornd(1,repmat(p_remain,1,p.Results.nsim));
remain   = remain > 0;

P.PP = P.PP(remain);
CP   = P;





