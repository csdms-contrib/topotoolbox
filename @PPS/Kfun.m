function [K,dvec] = Kfun(P,varargin)

%Kfun K-function of a point pattern on a stream network
%
% Syntax
%
%      [K,dvec] = Kfun(P)
%      [K,dvec] = Kfun(p,'pn',pv,...)
%
% Description
%
%      Kfun calculates the K-function of the point pattern on the stream
%      network P. Kfun uses the geometrically corrected K-function by Ang
%      et al. 2011.
%
% Input arguments
%
%      P     Instance of PPS
%
%      Parameter name/value pairs
%
%      'maxdist'   scalar
%      'val'   node attribute list with custom distance value. By default,
%              the distance is measured in mapunits
%      'usepoints' {false} or true: If false, the diameter is calculated
%              for the entire network, if true, only the distance between 
%              the points is calculated.
%      'plot'  {false} or true: If true, the function will create a plot
%              with the start and end point of the path with the maximum 
%              length.
%
% Output arguments
%
%      dmax    maximum distance 
%      IX      linear index into the GRIDobj with the end points of the
%              path with the maximum length.
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S,1);
%     P = PPS(S,'rpois',0.001,'z',DEM);
%     subplot(1,2,1)
%     [K,d] = Kfun(P,'maxdist',5000);
%     title('Random point pattern')
% 
%     P2 = generatepoints(P);
%     subplot(1,2,2)
%     [K,d] = Kfun(P2,'maxdist',5000);
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
addParameter(p,'nsim',100);
parse(p,varargin{:});

if ~isempty(p.Results.maxdist)
    maxdist = p.Results.maxdist;
else
    dmax = diameter(P,'val',p.Results.val,'d3d',p.Results.d3d);
    maxdist = dmax/20;
end

%% Simulate point patterns
total_length = tlength(P);
np = npoints(P);
C = simulate(P,'intensity',intensity(P) - 1/total_length,'nsim',p.Results.nsim);
IXSP = cellfun(@(x) x.PP,C,'UniformOutput',false);

%% calculate m
% For all points calculate the number of points that are exactly 
dvec = linspace(0,maxdist,p.Results.n+1);
ixgrid = P.S.IXgrid(P.PP);

Klocsim = cell(npoints(P),1);
Kloc = nan(npoints(P),numel(dvec)-1);

% Calculate local K functions
parfor r = 1:npoints(P)
    d = netdist(P.S,ixgrid(r));   
    % calculate m
    m = getlocation(P.S,dvec(2:end),'value',d,'output','cell');
    m = cellfun(@numel,m);

    % calculate actual local K function
    PP = P.PP(setdiff(1:npoints(P),r));
    n = histcounts(d(PP),dvec,'Normalization','cumcount');
    Kloc(r,:) = total_length/(np-1) * n./m;

    % local K function under CRS assumption
    nsim = nan(numel(C),p.Results.n);
    for rsim = 1:numel(IXSP)
        nsim(rsim,:) = histcounts(d(IXSP{rsim}),dvec(1:end),'Normalization','cumcount');
    end
    Klocsim{r} = total_length/(np-1) * nsim./m;
end

% the empirical K function
inan = isinf(Kloc) | isnan(Kloc);
Kloc(inan) = nan;
K = mean(Kloc,'omitnan');

% the simulated K function
Ksim = nan(p.Results.nsim,p.Results.n);
for r = 1:p.Results.nsim
    Kloc = cellfun(@(x) x(r,:),Klocsim,'UniformOutput',false);
    Kloc = vertcat(Kloc{:});
    inan = isinf(Kloc) | isnan(Kloc);
    Kloc(inan) = nan;
    Ksim(r,:) = mean(Kloc,'omitnan');
end
Ksim = [zeros(size(Ksim,1),1) Ksim];
Ksim = Ksim';

patch([dvec(:); flipud(dvec(:))],...
      [min(Ksim,[],2); flipud(max(Ksim,[],2))],...
      [.9 .9 .9],'EdgeColor','none');
hold on
plot([0 dvec(end)],[0 dvec(end)],'k--')
K = [0 K];
plot(dvec,K)
box on
xlim([0 dvec(end)])
hold off
xlabel('r')
ylabel('K(r)')
