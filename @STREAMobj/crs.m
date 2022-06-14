function [zs,exitflag,output] = crs(S,DEM,varargin)

%CRS constrained regularized smoothing of the channel length profile
%
% Syntax
%
%     zs = crs(S,DEM)
%     zs = crs(S,DEM,pn,pv,...)
%
% Description
%
%     Elevation values along stream networks are frequently affected by
%     large scatter, often as a result of data artifacts or errors. This
%     function returns a node attribute list of smoothed elevation values
%     calculated by nonparametric quantile regression with monotonicity
%     constraints. This function requires the Optimization Toolbox.
%
% Input parameters
%
%     S      STREAMobj
%     DEM    Digital elevation model (GRIDobj)
%     
%     Parameter name/value pairs {default}
%
%     'K'             positive scalar that dictates the degree of
%                     smoothing. The default is 2.
%     'tau'           quantile (scalar, >0 and <1). Default is 0.5.
%     'mingradient'   positive scalar {0}. Minimum downward gradient. 
%                     Choose carefully, because length profile may dip to
%                     steeply. Set this parameter to nan if you do not wish
%                     to have a monotonous dowstream decrease.
%     'nonstifftribs' {true} or false. true enables sharp knicks at
%                     confluences
%     'knickpoints'   nx1 vector with linear indices of locations
%                     where the stiffness penalty should be relaxed (e.g. 
%                     at knickpoints, waterfalls, weirs, dams, ...).
%                     The locations must be located on the stream network.
%                     This can be ensured by using STREAMobj/snap2stream.
%     'split'         {2}, 1 or 0. This option implements different ways to
%                     parallization. The option 2 successively processes
%                     trunk streams and may evaluate fastest. 1 processes
%                     individual drainage basins in parallel. 0 processes
%                     the entire network at one step. This may be the
%                     fastest option for small networks. Options 2 and 1
%                     may evaluate faster even without the availability of
%                     the parallel computing toolbox.
%
% Output parameters
%
%     zs     node attribute list with smoothed elevation values
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S);
%     S = trunk(S);
%     z10  = crs(S,DEM,'K',10,'tau',0.1);
%     z90  = crs(S,DEM,'K',10,'tau',0.9);
%     z50  = crs(S,DEM,'K',10,'tau',0.5);
%     plotdzshaded(S,[z10 z90])
%     hold on
%     plotdz(S,DEM,'color',[0.7 0.7 0.7])  
%     plotdz(S,z50,'color','k','LineWidth',2)
%     hold off
%     legend('K=10, tau=0.1-0.9 envelope',...
%            'original data',...
%            'K=10, tau=0.5')
%
% Algorithm
%
%     This algorithm uses nonparametric quantile regression to smooth the
%     data. The algorithm is described in Schwanghart and Scherler (2017)
%     (Eq. A13-A15).
%     
% References
%
%     Schwanghart, W., Scherler, D., 2017. Bumps in river profiles: 
%     uncertainty assessment and smoothing using quantile regression 
%     techniques. Earth Surface Dynamics, 5, 821-839. 
%     [DOI: 10.5194/esurf-5-821-2017]
%
% See also: STREAMobj/mincosthydrocon, quadprog, STREAMobj/crsapp, 
%           STREAMobj/crslin, STREAMobj/quantcarve, STREAMobj/smooth
%           
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 25. September, 2017


% check and parse inputs
narginchk(2,inf)

p = inputParser;
p.FunctionName = 'STREAMobj/crs';
addParamValue(p,'K',2,@(x) (isscalar(x) && x>0) || isa(x,'GRIDobj'));
addParamValue(p,'tau',0.5,@(x) isscalar(x) && x>0 && x<1)
addParamValue(p,'nonstifftribs',true,@(x) isscalar(x));
addParamValue(p,'mingradient',0,@(x) (isscalar(x) && (x>=0 || isnan(x))));
addParamValue(p,'knickpoints',[]);
addParamValue(p,'fixedoutlet',false);
addParamValue(p,'split',2);
parse(p,varargin{:});

% get node attribute list with elevation values
if isa(DEM,'GRIDobj')
    validatealignment(S,DEM);
    z = getnal(S,DEM);
elseif isnal(S,DEM)
    z = DEM;
else
    error('TopoToolbox:crs','Imcompatible format of second input argument')
end

if any(isnan(z))
    error('TopoToolbox:crs',...
        ['DEM or z may not contain any NaNs. Use STREAMobj/inpaintnans \n' ...
         'to fill nan values using interpolation'])
end

% double precision is necessary
z = double(z);

%% Run in parallel
if p.Results.split == 1
    % This option processes each connected component of the stream network
    % on a separate worker
    params = p.Results;
    params.split = false;
    [CS,locS] = STREAMobj2cell(S);
    if numel(CS) > 1
        % run only in parallel if more than one drainage basin
        Cz = cellfun(@(ix) z(ix),locS,'UniformOutput',false);
        Czs = cell(size(CS));
        parfor r = 1:numel(CS)
            Czs{r} = crs(CS{r},Cz{r},params);
        end
        
        zs = nan(size(z));
        for r = 1:numel(CS)
            zs(locS{r}) = Czs{r};
        end
        return
    else
        zs = crs(S,z,params);
        return
    end
    
elseif p.Results.split == 2
    % This option processes each tributary on a separate worker.
    % Tributaries are recursively extracted from the stream network,
    % starting with the trunk river, then with the longest rivers tributary
    % to the trunk river, and so on. 
        
    [CS,locb,CID] = STREAMobj2cell(S,'trib');
    
    params = p.Results;
    params.split = 0;
    % params.fixedoutlet = false;
    for r = 1:max(CID)
        ii = CID == r;
        if r > 1
            params.fixedoutlet = true;
        end
        CStemp   = CS(ii);
        locbtemp = locb(ii);
        ztribs    = cellfun(@(ix) z(ix),locbtemp,'UniformOutput',false);
        Czstemp   = cell(numel(ii),1);
        
        parfor r2 = 1:numel(CStemp)
            Czstemp{r2} = crs(CStemp{r2},ztribs{r2},params);
        end
        for r2 = 1:numel(CStemp)
            z(locbtemp{r2}) = Czstemp{r2};
        end
    end
    
    zs = z;
    return
        
    
end

%% CRS starts here
% upstream distance
d  = S.distance;

% nr of nodes
nr = numel(S.IXgrid);

if nr < 3
    % a special case. The stream network needs at least three nodes in a row. If
    % there are less, the second derivative matrix is empty, and crs uses
    % quantile carving instead
    zs = quantcarve(S,z,p.Results.tau,...
                        'split',p.Results.split,...
                        'mingradient',p.Results.mingradient,...
                        'fixedoutlet',p.Results.fixedoutlet);
    return
end


%% Second-derivate matrix
%
% This matrix contains the finite central differences 
[I,loc] = ismember(S.ixc,S.ix);

%         i-1           i        i+1
colix  = [S.ixc(loc(I)) S.ixc(I) S.ix(I)];

% Set user-supplied knickpoints to non-stiff
if ~isempty(p.Results.knickpoints)
    IX = p.Results.knickpoints;
    I  = ismember(S.IXgrid,IX);
    I  = I(colix(:,2));
    colix(I,:) = [];
end

val    = [2./((d(colix(:,2))-d(colix(:,1))).*(d(colix(:,3))-d(colix(:,1)))) ...
         -2./((d(colix(:,3))-d(colix(:,2))).*(d(colix(:,2))-d(colix(:,1)))) ...
          2./((d(colix(:,3))-d(colix(:,2))).*(d(colix(:,3))-d(colix(:,1))))];

% Set tributaries to non-stiff
if p.Results.nonstifftribs
    dd = distance(S,'max_from_ch');
    I  = (dd(colix(:,2)) - dd(colix(:,3)))>=(sqrt(2*S.cellsize.^2)+S.cellsize/2);
    colix(I,:) = [];
    val(I,:) = [];
end

% second-derivative matrix
nrrows = size(colix,1);
rowix  = repmat((1:nrrows)',1,3);
Asd    = sparse(rowix(:),colix(:),val(:),nrrows,nr);

if isempty(Asd)
    % a special case. The stream network needs at least three nodes in a row. If
    % there are less, the second derivative matrix is empty, and crs uses
    % quantile carving instead
    zs = quantcarve(S,z,p.Results.tau,...
                        'split',p.Results.split,...
                        'mingradient',p.Results.mingradient,...
                        'fixedoutlet',p.Results.fixedoutlet);
    return
end


%% Setup linear system
% balance stiffness and fidelity
Asd    = p.Results.K * S.cellsize^2 * sqrt(nr/nrrows) * Asd;
C      = [sparse(nr,nr);Asd];

%% Monotonicity constraint
if ~isnan(p.Results.mingradient)
    d = 1./(d(S.ix)-d(S.ixc));
    A = [sparse(nr,nr*2) (sparse(S.ix,S.ixc,d,nr,nr)-sparse(S.ix,S.ix,d,nr,nr))];
    e = zeros(nr,1);
    if p.Results.mingradient > 0
        e(S.ix) = -p.Results.mingradient;
    end
else
    A = [];
    e = [];
end

%% reformulate for quadprog
H = 2*(C'*C);
H = [sparse(nr*2,nr*3);sparse(nr,nr*2) H]; 

b = [z; zeros(nrrows,1)];
c = -2*C'*b;
f = [p.Results.tau*ones(nr,1);(1-p.Results.tau)*ones(nr,1);c];

% Equalities
if ~p.Results.fixedoutlet
    Aeq = [speye(nr),-speye(nr),speye(nr)];
else 
    OUTL = streampoi(S,'outlet','logical');
    P    = spdiags(+(~OUTL),0,nr,nr);
    Aeq  = [P,-P,speye(nr)];
end
beq = z;

% Lower and upper bounds
lb  = [zeros(nr,1);zeros(nr,1);-inf*ones(nr,1)];
ub  = inf(nr*3,1);
    
% options
if verLessThan('optim','6.3')
    options = optimset('Display','off','Algorithm','interior-point-convex');
else
    options = optimoptions('quadprog','Display','off');
end

%% Solve
[zs,~,exitflag,output] = quadprog(H,f,A,e,Aeq,beq,lb,ub,[],options);

zs = zs(2*nr+1:end);

end

