function OUT = Kfun(P,varargin)

%Kfun K-function of a point pattern on a stream network
%
% Syntax
%
%      K = Kfun(P)
%      K = Kfun(p,'pn',pv,...)
%
% Description
%
%      Kfun calculates the K-function of the point pattern on the stream
%      network P. Kfun uses the geometrically corrected K-function by Ang
%      et al. 2011, so that the theoretical K-function of a completely
%      random point process is K(r) = r. For large distances, however, the 
%      geometric correction, which is similar to Ripleys edge correction
%      of a 2D point pattern, may not work properly.
%
% Input arguments
%
%      P     Instance of PPS
%
%      Parameter name/value pairs
%
%      'method' Either 'ang' or 'okabe'. 'ang' is the geometrically
%              corrected K function while 'okabe' is uncorrected. This 
%              means that there is no theoretic CRS K function for 'okabe'.
%      'maxdist'   scalar. Maximum distance used to calculate the K-
%              function. The default is diameter(P)/20.
%      'val'   node attribute list with custom distance value. By default,
%              the distance is measured in mapunits
%      'plot'  {false} or true: If true, the function will create a plot
%              with the start and end point of the path with the maximum 
%              length.
%      'int'   node attribute list with intensity values. If provided then
%              Kfun calculates the inhomogeneous K-function (currently, 
%              this is only supported for method 'okabe'. 
%
% Output arguments
%
%      K      structure array with following fields
%      .d     distance vector
%      .K     empirical K values
%      .Ktheo theoretical K values (if method is 'ang')
%      .KsimU Upper bound on acceptance intervals of the simulated K
%             function
%      .KsimL Lower bound
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S,1);
%     P = PPS(S,'rpois',0.001,'z',DEM);
%     subplot(1,2,1)
%     K1 = Kfun(P,'maxdist',5000);
%     title('Random point pattern')
% 
%     P2 = generatepoints(P);
%     subplot(1,2,2)
%     K2 = Kfun(P2,'maxdist',5000);
%     title('Regular point pattern')
%
% See also: PPS
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 22. December, 2020


%% Calculate graph diameter to get maximum distance

p  = inputParser;
addParameter(p,'d3d',false);
addParameter(p,'val',[]);
addParameter(p,'maxdist',[]);
addParameter(p,'n',50);
addParameter(p,'nsim',39);
addParameter(p,'method','ang');
addParameter(p,'plot',true);
addParameter(p,'int',[],@(x) isnal(P.S,x));
parse(p,varargin{:});

if ~isempty(p.Results.maxdist)
    maxdist = p.Results.maxdist;
else
    dmax = diameter(P,'val',p.Results.val,'d3d',p.Results.d3d);
    maxdist = dmax/20;
end

nsim = p.Results.nsim;
method = validatestring(p.Results.method,{'ang','okabe'});

%% Simulate point patterns
total_length = tlength(P);
np = npoints(P);

if nsim > 0
    if isempty(p.Results.int)
        C = simulate(P,'intensity',intensity(P) - 1/total_length,'nsim',p.Results.nsim);
    else
        int = p.Results.int * (1-1/total_length);
        C = simulate(P,'intensity',int,'nsim',p.Results.nsim);
    end
    IXSP = cellfun(@(x) x.PP,C,'UniformOutput',false);
end

%% calculate m
% For all points calculate the number of points that are exactly 
dvec = linspace(0,maxdist,p.Results.n+1);
ixgrid = P.S.IXgrid(P.PP);

if nsim > 0
Klocsim = cell(npoints(P),1);
else
    C    = [];
    IXSP = [];
end
Kloc = nan(npoints(P),numel(dvec)-1);

% Calculate local K functions

if isempty(p.Results.int)
    intp = [];
    int  = [];
else
    intp = getmarks(P,p.Results.int);
    int  = p.Results.int;
end
    

for r = 1:npoints(P)
    d = netdist(P.S,ixgrid(r));
    m = [];
    % calculate m
    switch method
        case 'ang'
            % get number of locations that are exactly a distance r away
            % from the point
            m = getlocation(P.S,dvec(2:end),'value',d,'output','cell');
            m = cellfun(@numel,m);
    end

    % calculate actual local K function
    PP = P.PP(setdiff(1:npoints(P),r));
    
    switch method
        case 'ang'
            if isempty(intp)
                n = histcounts(d(PP),dvec,'Normalization','cumcount');
                Kloc(r,:) = total_length/(np-1) * n./m;
            else
                [~,~,bin] = histcounts(d(PP),dvec,'Normalization','count');
                n = accumarray(bin,intp(PP).*m(PP),[],@(x) sum(1./(intp(r)*x)));
                n = cumsum(n);
                Kloc(r,:) = n;
            end
                
        case 'okabe'
            n = histcounts(d(PP),dvec,'Normalization','cumcount');
            Kloc(r,:) = n;
    end

    if nsim > 0
        % local K function under CRS assumption
        sim = nan(numel(C),p.Results.n);
        if isempty(intp) || strcmp(method, 'okabe')
            for rsim = 1:numel(IXSP)
                sim(rsim,:) = histcounts(d(IXSP{rsim}),dvec(1:end),'Normalization','cumcount');
            end
        else
            for rsim = 1:numel(IXSP)
                [~,~,bin] = histcounts(d(PP),dvec,'Normalization','count');
                intsim = int(IXSP{rsim});
                n = accumarray(bin,intsim(PP).*m(PP),[],@(x) sum(1./(intp(r)*x)));
                n = cumsum(n);
                sim(rsim,:) = n;
            end
        end
                
        
        switch method
            case 'ang'
                if isempty(intp)
                    Klocsim{r} = total_length/(np-1) * sim./m;
                else
                    Klocsim{r} = sim;
                end
                    
            case 'okabe'
                Klocsim{r} = sim;
        end
    end
end

% the empirical K function
inan = isinf(Kloc) | isnan(Kloc);
Kloc(inan) = nan;
switch method
    case 'ang'
        if isempty(p.Results.int)
            K = mean(Kloc,'omitnan');
        else
            K = 1/sum(1./intp) * sum(Kloc,'omitnan');
        end
    case 'okabe'
        K = sum(Kloc,'omitnan');
        K = total_length/(np*(np-1)).*K;
end

% the simulated K function
if nsim > 0
Ksim = nan(p.Results.nsim,p.Results.n);
for r = 1:p.Results.nsim
    Kloc = cellfun(@(x) x(r,:),Klocsim,'UniformOutput',false);
    Kloc = vertcat(Kloc{:});
    inan = isinf(Kloc) | isnan(Kloc);
    Kloc(inan) = nan;
    switch method
        case 'ang'
            if isempty(p.Results.int)
                Ksim(r,:) = mean(Kloc,'omitnan');
            else
                Ksim(r,:) = 1/sum(1./intp) * sum(Kloc,'omitnan');
            end
        case 'okabe'
            Ksim(r,:) = sum(Kloc,'omitnan');
            Ksim(r,:) = total_length/(np*(np-1)).*Ksim(r,:);
    end
end

Ksim = [zeros(size(Ksim,1),1) Ksim];
Ksim = Ksim';

patch([dvec(:); flipud(dvec(:))],...
      [min(Ksim,[],2); flipud(max(Ksim,[],2))],...
      [.9 .9 .9],'EdgeColor','none');
hold on
end

if p.Results.plot
switch method
    case 'ang'
        plot([0 dvec(end)],[0 dvec(end)],'k--')
    case 'okabe'
        plot(dvec,mean(Ksim,2),'k--')
end
hold on
K = [0 K];
plot(dvec,K)
box on

xlim([0 dvec(end)])
hold off
xlabel('r [m]')
switch method
    case 'ang'
        ylabel('K_l(r)')
    case 'okabe'
        ylabel('K(r)')
end
end
if nargout >= 1
    OUT.d = dvec;
    OUT.K = K;
    
    switch method
        case 'ang'
            OUT.Ktheo = dvec;
        case 'okabe'
            OUT.Kmeansim = mean(Ksim,2);
    end
    
    if nsim > 0
        OUT.KsimU = max(Ksim,[],2);
        OUT.KsimL = min(Ksim,[],2);
    end
end
    
