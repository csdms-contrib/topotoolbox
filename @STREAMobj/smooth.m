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
%     Smooth can handle nan values and takes weights as additional input
%     arguments. 
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
%                the parallel processing toolbox). If there is only one
%                drainage basin, there is no gain in setting split to true.
%
%     Only applicable if method = 'regularization'
%
%     'K'        scalar, stiffness penalty (positive scalar, default = 10)
%     'nstribs'  relax stiffness penalty at tributary junctions (true
%                (default) or false)
%     'positive' set true if zs must be positive (uses lsqlin and requires 
%                the optimization toolbox)
%     'weights'  GRIDobj or node-attribute list with weights. Weights must 
%                be >= 0. By default, all weights are one. 
%     'breaks'   vector with linear indices into the GRIDobj from which S
%                was derived, or PPS object. The smoothness constraint is
%                relaxed at these points.
%
% Output parameters
%
%     zs     node attribute list with smoothed elevation values
%
% Example 1
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
% Example 2
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',100);
%     S = klargestconncomps(S);   
%     St = trunk(S);
%     G  = gradient8(DEM);
%
%     % Upstream hillslope area
%     A  = upslopestats(FD,GRIDobj(DEM)+1,'sum',S);
%     G  = upslopestats(FD,G,'mean',S);
%     
%     % Derive and plot mean hillslope angle along 
%     % the trunk river
%     plotdz(St,DEM);
%     yyaxis right
%     scatter(St.distance,getnal(St,G),getnal(St,A),'ko')
%     gs = smooth(St,G,'K',1000);
%     hold on
%     plotdz(St,gs)
%     gsw = smooth(St,G,'K',1000,'weights',A);
%     plotdz(St,gsw,'LineWidth',2)
%     ylabel('Gradient [-]')
%     hold off
%     legend('River profile',...
%           'Hillslope gradient',...
%           'Smoothed hillslope gradient',...
%           'Weighted smoothed hillslope gradient')
%
% Algorithm
%
%     This algorithm uses regularized interpolation to smooth the data. The
%     algorithm is described in Schwanghart and Scherler (2017) (Eq. A6-A10). 
%     Setting the 'positive' to true forces the returned values to be
%     positive. In this case, the function uses nsqlin with a lower bound
%     of zero for all nodes.
%
% References
%
%     Schwanghart, W., Scherler, D., 2017. Bumps in river profiles: 
%     uncertainty assessment and smoothing using quantile regression 
%     techniques. Earth Surface Dynamics, 5, 821-839. 
%     [DOI: 10.5194/esurf-5-821-2017]
%
% See also: STREAMobj/crs, STREAMobj/crsapp, STREAMobj/inpaintnans,
%           STREAMobj/quantcarve, STREAMobj/crslin
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 2. September, 2020


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
addParameter(p,'weights',[],@(x) isa(x,'GRIDobj') || isnal(S,x) || isempty(x));
addParameter(p,'breaks',[]);
addParameter(p,'distance',[]);

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

% handle weights
weights = p.Results.weights;
if isempty(weights) 
elseif isa(weights,'GRIDobj')
    validatealignment(S,weights);
    weights = getnal(S,weights);
elseif isnal(S,weights)
    % great, go on
end

% handle custom distances
if isnal(S,p.Results.distance)
    d = p.Results.distance;
end

% check for nans
if strcmp(p.Results.method,'regularization') && any(isnan(z))
    % error('DEM or z may not contain any NaNs if method is regularization.')
    inan = isnan(z);
    z = inpaintnans(S,z,'extrap',true);
    if isempty(weights)
        weights = double(~inan);   
    else
        weights(inan) = 0;
    end
end

% run in parallel if wanted
if p.Results.split
    params = p.Results;
    params.split = false;
    [CS,locS] = STREAMobj2cell(S);
    Cz = cellfun(@(ix) z(ix),locS,'UniformOutput',false);
    Czs = cell(size(CS));
    
    % handle weights
    if ~isempty(params.weights)
        Cw = cellfun(@(ix) weights(ix),locS,'UniformOutput',false);
    else
        Cw = cell(size(CS));
    end
    params = rmfield(params,'weights');
    
    % handle custom distances
    if ~isempty(params.distance)
        Cd = cellfun(@(ix) d(ix),locS,'UniformOutput',false);
    else
        Cd = cell(size(CS));
    end
    params = rmfield(params,'distance');
    
    % finally, do the smoothing
    parfor r = 1:numel(CS)
        Czs{r} = smooth(CS{r},Cz{r},params,'weights',Cw{r},'distance',Cd{r});
    end

    % write smoothed values to node-attribute list
    zs = nan(size(z));
    for r = 1:numel(CS)
        zs(locS{r}) = Czs{r};
    end
    return
end

%% Smoothing starts here
% upstream distance
if ~isempty(p.Results.distance)
    d = p.Results.distance;
else
    d = S.distance;
end

% d  = S.distance;
% nr of nodes
nr = numel(S.IXgrid);

if nr <= 2
    % If there are no more than 2 nodes in the network, we won't smooth
    zs = z;
    return
end

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
        
        if ~isempty(p.Results.breaks)
            % do breaks come as a PPS object
            if isa(p.Results.breaks,'PPS')
                kp = points(p.Results.breaks,'IXgrid');
            else
                kp = p.Results.breaks;
            end
            
            % identify breaks that belong to the basin 
            I  = ismember(S.IXgrid,kp);
            
            % remove stiffness at break nodes
            I  = I(colix(:,2));
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
        b      = double(b);
        
        % use weighted least squares if required
        % weights are normalized so that their sum equals the number
        % of rows in the fidelity matrix. 
        if ~isempty(weights)
            weights = double(weights);
            weights = weights./sum(weights) * nr;
            W = spdiags([weights;ones(nrrows,1)],0,nr+nrrows,nr+nrrows);
            C = W*C;
            b = W*b;
        end
        
        if ~p.Results.positive
            zs     = C\b;
        else
            options = optimset('Display','off');
            zs = lsqlin(C,b,[],[],[],[],zeros(nr,1),[],[],options);

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
        
