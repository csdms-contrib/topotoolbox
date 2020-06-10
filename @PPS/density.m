function [h,bw] = density(P,varargin)

%DENSITY nonparametric density estimation on network
%
% Syntax
%
%     h = density(P)
%     h = density(P,'bw',bandwidth)
%     h = density(P,'bw',[],...)
%     
% Description
%
%     density calculates a density estimation of the points in point
%     pattern P. It is a measure of how meany points are found in the
%     vicinity of each location in the stream network.
%
%     The function uses a diffusion estimator (McSwiggan et al., 2017)
%     which corresponds to an equal-split continuous Gaussian Kernel. Other
%     than McSwiggan et al., the function uses an implicit solver of the
%     diffusion equation.
%     
%     The bandwidth is supplied as optional parameter name/value 'bw'. If
%     this argument is empty the function attempts to find an optimal
%     bandwidth using cross-validation. This, however, can take some time
%     and may sometimes fail. Bandwidth optimization is done using Bayesian
%     Optimization.
%
% Input arguments
%
%     P      point pattern on stream network (class PPS)
%
%     Parameter name/value pairs
%
%     'bandwidth' bandwidth (scalar) in map units. The default is 20*cellsize.
%            If you set the parameter 'distance', then 'bandwidth' must be in the
%            units of distance. If the bandwidth is supplied as empty
%            array, the function will compute a cross-validated estimate of
%            an optimal bandwidth using Bayesian optimization. This
%            computation can be quite long.
%     'nriter'  Though the implicit solution is unconditionally stable,
%            the precision of the heat kernel can be improved by a larger 
%            number of iterations. The default is 10.
%     'cvsamples' Number of points used for cross-validation. Decreasing 
%            the number of samples decreases runtime of the optimization at
%            the cost of a random objective function.
%     'bwrange' two-element vector with minimum and maximum values for the
%            bandwidth cross-validation. 
%     'useparallel' {true} or false.
%     'weights' weights for each point in P. Default is ones(npoints(P),1).
%            Only positive values or zeros.
%     'distance' node attribute list of distances (default is P.S.distance)
%     'bc'   Boundary conditions for the diffusion estimator. 
%            0 = all open boundaries (default)
%            1 = outlet open only 
%            2 = channelheads open only
%
% Output arguments
%
%     h     density (node-attribute list)
%     bw    bandwidth
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','c');
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S,1);
%     rng(2); % make results replicable
%     P = PPS(S,'runif',3,'z',DEM);
%     d = density(P,'bw',2000);
%     P = simulate(P,'intensity',d*15);
%     [d,bw] = density(P,'bw',[],'bwrange',[100 4000]);
%     plotc(S,d)
%     hold on
%     plotpoints(P)
%     hold off
%
% Reference
%
%     McSwiggan G, Baddeley A, Nair G. 2017. Kernel Density Estimation on a
%     Linear Network. Scandinavian Journal of Statistics 44 : 324–345. DOI:
%     10.1111/sjos.12255
%
% See also: PPS, PPS/npoints, PPS/intensity, PPS/histogram, PPS/rhohat
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 23. September, 2019

p = inputParser;
p.FunctionName = 'PPS/density';
addParameter(p,'bandwidth',20*P.S.cellsize)
addParameter(p,'nriter',10);
addParameter(p,'cvsamples',100);
addParameter(p,'useparallel',true);
addParameter(p,'bwrange',[P.S.cellsize,P.S.cellsize*100]);
addParameter(p,'distance',[],@(x) isnal(P.S,x));
addParameter(p,'weights',ones(size(P.PP)),@(x) (numel(x) == npoints(P)) && all(x>=0));
addParameter(p,'bc',0);
parse(p,varargin{:});

% If there are no points
if npoints(P) == 0
    bw = [];
    h  = getnal(P.S);
    return
end

S = P.S;
% Edge distances
if isempty(p.Results.distance)
    d = sqrt((S.x(S.ix) - S.x(S.ixc)).^2 + (S.y(S.ix) - S.y(S.ixc)).^2);
else
    if ~isnal(S,p.Results.distance)
        error('Distance is not a node-attribute list.');
    end
    d = abs(p.Results.distance(S.ix) - p.Results.distance(S.ixc));
end

% Cross-validation samples
cvsamples = min(p.Results.cvsamples,npoints(P));

% Boundary conditions
if p.Results.bc == 0
    endnodes  = streampoi(S,{'chan','outlet'},'logical');
elseif p.Results.bc == 1
    endnodes  = streampoi(S,'outlet','logical');
elseif p.Results.bc == 2
    endnodes  = streampoi(S,'chan','logical');
else
    error('Wrong value for bc')
end

% Adjacency matrix
n = numel(S.x);
A = sparse(S.ix,S.ixc,1./(d.^2),n,n);
A = A+A';

% Points with weights

nnodes = numel(P.S.x);
w      = p.Results.weights./sum(p.Results.weights) .* npoints(P);
h      = accumarray(P.PP,w(:),[nnodes 1],@sum,0);

if ~isempty(p.Results.bandwidth)
    bw     = p.Results.bandwidth;
    t      = bw.^2;
    nriter = p.Results.nriter;
    dt     = t/nriter;
    A      = 1/2 * dt * A;
else
    vars = optimizableVariable('bw',p.Results.bwrange,'Type','real');
    t      = bayesopt(@(x) crossvalD(x),vars,'UseParallel',p.Results.useparallel);
    bP     = bestPoint(t);
    bw     = bP.bw;
    t      = bw.^2;
    nriter = p.Results.nriter;
    dt     = t/nriter;
    A      = 1/2 * dt * A;
    
end

% Laplacian
L = -A + spdiags(endnodes+1,0,n,n) + spdiags(sum(A,1)',0,n,n);

% L = -A + speye(n) + alpha*spdiags(sum(A,2),0,n,n);

% implicit
h = diffusion_impl(L,h,nriter);
h = h./P.S.cellsize;
% h = h./(P.S.cellsize.^2);


function negll = crossvalD(vars)
    
    bw = vars.bw;
    
    t      = bw.^2;
    nriter = p.Results.nriter;
    dt     = t/nriter;
    B      = 1/2 * dt * A;
    
    L     = -B + spdiags(endnodes+1,0,n,n) + spdiags(sum(B,1)',0,n,n);
    ht    = h;
    ht    = diffusion_impl(L,ht,nriter);
    ppp   = P.PP(randperm(npoints(P)));
    ll    = zeros(cvsamples,1);
    parfor rr = 1:cvsamples
        pc = 1:numel(ppp);
        pc(rr) = [];
        % pc = setdiff(1:numel(ppp),rr);
        h2 = zeros(n,1);
        h2(ppp(pc)) = 1;

        h2 = diffusion_impl(L,h2,nriter);
        ll(rr) = h2(ppp(rr));
    end

    negll = sum((ll-ht(ppp(1:cvsamples))).^2);
end
    
end

function h = diffusion_impl(L,h,nriter)
for iter = 1:nriter
    h = L'\h;
end
end



