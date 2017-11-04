function zs = smooth(S,DEM,varargin)

%SMOOTH smoothing of node-attribute lists
%
% Syntax
%
%     zs = smooth(S,DEM)
%     zs = smooth(S,z)
%     zs = smooth(...,pn,pv,...)
%
% Description
%
%     Elevation values along stream networks are frequently affected by
%     large scatter, often as a result of data artifacts or errors. This
%     function returns a node attribute list of elevations calculated by
%     smoothing. For further control on the results use the function
%     STREAMobj/crs.
%
% Input parameters
%
%     S        STREAMobj
%     DEM      Digital elevation model (GRIDobj) or any other GRIDobj
%     z        node attribute list
%
% Parameter name/value pairs
%
%     'method'   {'regularization'} or 'movmean'. 'movmean' calculates a
%                three-point moving average.
%     'split'    {false} or true. True will identify individual drainage
%                basins and process each individually in parallel (requires
%                the parallel processing toolbox).
%
%     Only applicable if method = 'regularization'
%
%     'K'        scalar, stiffness penalty (positive scalar, default = 10)
%     'nstribs'  relax stiffness penalty at tributary junctions (true
%                (default) or false)
%     'positive' set true if zs must be positive (uses lsqlin and requires 
%                the optimization toolbox)
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
%     g = gradient(S,imposemin(S,DEM));
%     gs = smooth(S,g,'k',100);
%     subplot(2,1,1)
%     plotc(S,g)
%     box on
%     colorbar
%     axis image
%     subplot(2,1,2)
%     plotc(S,gs)
%     axis image
%     box on
%     colorbar
%
% Algorithm
%
%     This algorithm uses regularized interpolation to smooth the data.
%     Setting the 'positive' to true forces the returned values to be
%     positive. In this case, the function uses nsqlin with a lower bound
%     of zero for all nodes.
%
% See also: STREAMobj/crs, STREAMobj/crsapp, STREAMobj/inpaintnans,
%           STREAMobj/quantcarve, STREAMobj/crslin
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 18. July, 2017


% check and parse inputs
narginchk(2,inf)

p = inputParser;
p.FunctionName = 'STREAMobj/smooth';
addParameter(p,'method','regularization'); 
addParameter(p,'split',false);
% parameters for regularization
addParameter(p,'K',10,@(x) (isscalar(x) && x>0));
addParameter(p,'nstribs',true,@(x) isscalar(x));
addParameter(p,'positive',false,@(x) isscalar(x));

parse(p,varargin{:});

method = validatestring(p.Results.method,{'regularization','movmean'});

% get node attribute list with elevation values
if isa(DEM,'GRIDobj')
    validatealignment(S,DEM);
    z = getnal(S,DEM);
elseif isnal(S,DEM)
    z = DEM;
else
    error('Imcompatible format of second input argument')
end

% check for nans
if strcmp(p.Results.method,'regularization') && any(isnan(z))
    error('DEM or z may not contain any NaNs if method is regularization.')
end

% run in parallel if wanted
if p.Results.split
    params = p.Results;
    params.split = false;
    [CS,locS] = STREAMobj2cell(S);
    Cz = cellfun(@(ix) z(ix),locS,'UniformOutput',false);
    Czs = cell(size(CS));
    parfor r = 1:numel(CS)
        Czs{r} = smooth(CS{r},Cz{r},params);
    end
    
    zs = nan(size(z));
    for r = 1:numel(CS)
        zs(locS{r}) = Czs{r};
    end
    return
end

%% Smoothing starts here
% upstream distance
d  = S.distance;
% nr of nodes
nr = numel(S.IXgrid);

switch method
    case 'regularization'

        %% Fidelity matrix
        %
        % This matrix is an identity matrix, since we are only interested in 
        % predicting values at the node locations of the network.
        Afid = speye(nr,nr);      

        %% Second-derivate matrix
        %
        % This matrix contains the finite central differences 
        [I,loc] = ismember(S.ixc,S.ix);

        %         i-1           i        i+1
        colix  = [S.ixc(loc(I)) S.ixc(I) S.ix(I)];

        val    = [2./((d(colix(:,2))-d(colix(:,1))).*(d(colix(:,3))-d(colix(:,1)))) ...
                 -2./((d(colix(:,3))-d(colix(:,2))).*(d(colix(:,2))-d(colix(:,1)))) ...
                  2./((d(colix(:,3))-d(colix(:,2))).*(d(colix(:,3))-d(colix(:,1))))];

        % Set tributaries to non-stiff
        if p.Results.nstribs
            dd = distance(S,'max_from_ch');
            I = (dd(colix(:,2)) - dd(colix(:,3)))>=(sqrt(2*S.cellsize.^2)+S.cellsize/2);
            colix(I,:) = [];
            val(I,:) = [];
        end

        % second-derivative matrix
        nrrows = size(colix,1);
        rowix  = repmat((1:nrrows)',1,3);
        Asd    = sparse(rowix(:),colix(:),val(:),nrrows,nr);

        %% Setup linear system
        % balance stiffness and fidelity
        F      = p.Results.K * sqrt(size(Afid,1)/nrrows) * S.cellsize^2;
        % F      = p.Results.K * norm(Afid,1)/norm(Asd,1);
        Asd    = F*Asd;
        C      = [Afid;Asd];
        % right hand side of equation
        b      = [z;zeros(nrrows,1)]; % convex -> negative values

        % no minimum gradient imposition involves only solving the
        % overdetermined system of equations. 
        if ~p.Results.positive
            zs     = C\b;
        else
            zs = lsqlin(C,b,[],[],[],[],zeros(nr,1),[],z);
        end
        
    case 'movmean'
        
        % movmean implements a moving three-point average. It calculates
        % the upstream and downstream gradients from each point and
        % averages the upstream gradients since there may be more than one
        % upstream neighbor. 
        
        g  = gradient(S,z);
        cs = S.cellsize;
        zs = (3*z - g*cs + accumarray(S.ixc,g(S.ix),size(z),@(x) mean(x)*cs))/3;

end

end
        
