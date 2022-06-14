function varargout = ecdf(P,varargin)

%ECDF Empirical cumulative distribution function of a covariate 
%
% Syntax
%
%     ecdf(P)
%     ecdf(P,pn,pv,...)
%     [CDF,x,CDFpoi,CDFlower,CDFupper,bw] = ...
% 
% Description
%
%     ecdf returns the empirical cumulative distribution function (eCDF) of the
%     point pattern P. Without further input arguments, ecdf computes the 
%     stepwise function based on the stream network distance calculated
%     from the outlets of the stream network in P.
%
%     ecdf compares the computed eCDF to the hypothesis of complete spatial
%     randomness (CSR) and calculates either acceptance or confidence
%     intervals using bootstrapping.
%
% Input arguments
%
%     P      point pattern on stream network (class PPS)
%
%     Parameter name/value pairs
%
%     'covariate'     Node-attribute list of covariate values. By default, 
%                     the distance from the outlet is used. Alternatives are
%                     'z' or e.g. gradient(P.S,DEM).
%     'accintervals'  {true} or false. true plots acceptance intervals.
%                     Acceptance intervals are the bounds around the
%                     hypothesized CSR model (complete spatial randomness). 
%                     The width reflects the inherent uncertainty of the
%                     data if the hypothesis of the model is true.
%     'confintervals' true or {false}. true plots confidence intervals.
%                     Confidence intervals are bounds around the eCDF
%                     curve and reflects the confidence around the
%                     calculated eCDF. 
%     'nsim'          Number of bootstrap simulations (1000)
%     'plot'          {true} or false
%     'alpha'         width of confidence intervals (1-alpha). Default is 
%                     0.05.
%     'ksdensity'     true or {false}. Uses kernel density estimates using
%                     the function ksdensity. 
%     'bandwidth'     bandwidth of kernel density. Only applicable if
%                     'ksdensity' is true.
%
% Output arguments
%
%     CDF         values of eCDF
%     x           values on x-axis of eCDF
%     CDFpoi      hypothesized CDF (CSR)
%     CDFlower    lower and
%     CDFupper    upper bounds of the confidence or acceptance limits
%     bw          bandwidth selected if 'ksdensity' is true
%
% Example
%
%      DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%      FD  = FLOWobj(DEM,'preprocess','c');
%      S = STREAMobj(FD,'minarea',1000);
%      S = klargestconncomps(S,1);
%      P = PPS(S,'runif',100,'z',DEM);
%      subplot(1,2,1)
%      ecdf(P,'accintervals',true)
%      subplot(1,2,2)
%      ecdf(P,'confintervals',true)
%
% See also: PPS, PPS/rhohat
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. February, 2019



% Check input arguments
p = inputParser;
p.FunctionName = 'ecdf';
addRequired(p,'P',@(x) isa(x,'PPS'));

% Add elevation
addParameter(p,'accintervals',true);
addParameter(p,'confintervals',false);
addParameter(p,'nsim',1000);
addParameter(p,'covariate',[]);
addParameter(p,'plot',true);
addParameter(p,'alpha',0.05);
addParameter(p,'ksdensity',false);
addParameter(p,'bandwidth',[]);
addParameter(p,'function','cdf');
addParameter(p,'weights',ones(npoints(P),1),@(x) numel(x) == npoints(P) & all(x > 0));

% Parse
parse(p,P,varargin{:});

% get covariate
c = getcovariate(P,p.Results.covariate);

% x-axis
edges = linspace(min(c),max(c),1000);
x     = edges(1:end-1)+diff(edges([1 2]));


useks  = p.Results.ksdensity;
usefun = p.Results.function;

if ~useks
    % histogram
    N  = histcounts(c(P.PP),edges,'Normalization',usefun);
    Np = histcounts(c,edges,'Normalization',usefun);
    bw = 0;
else
    % kernel density
    [N,~,bw] = ksdensity(c(P.PP),x,'function',usefun,'bandwidth',p.Results.bandwidth,'weights',p.Results.weights);
    Np = ksdensity(c,x,'function',usefun,'bandwidth',bw);
end

nnp   = numel(P.PP);
if p.Results.accintervals && ~p.Results.confintervals
    NS    = zeros(numel(edges)-1,p.Results.nsim);
    S     = P.S;
    parfor r = 1:p.Results.nsim
        PS = PPS(S,'runif',nnp);
        csim = PS.PP;
        if ~useks
            NS(:,r) = histcounts(c(csim),edges,'Normalization',usefun);
        else
            NS(:,r) = ksdensity(c(csim),x,'function',usefun,'bandwidth',bw);
        end
    end
    legendentry = 'Acceptance interval';
end

if p.Results.confintervals
    % create bootstrap samples
    RI   = randi(nnp,nnp,p.Results.nsim);
    csim = c(P.PP(RI));
    wsim = p.Results.weights(RI);
    NS    = zeros(numel(edges)-1,p.Results.nsim);
    parfor r = 1:p.Results.nsim
        if ~useks
            NS(:,r)   = histcounts(csim(:,r),edges,'Normalization',usefun);
        else
            NS(:,r)   = ksdensity(csim(:,r),x,'function',usefun,...
                'bandwidth',bw,'weights',wsim(:,r));
        end
    end    
    legendentry = 'Confidence interval';
end

% Calculate bounds
if p.Results.accintervals || p.Results.confintervals
NQ   = quantile(NS,[p.Results.alpha 1-p.Results.alpha],2);
clear NS
end

x = x(:);

if p.Results.plot
    tf = ishold;
    if p.Results.accintervals || p.Results.confintervals
        patch([x; flipud(x)],[NQ(:,1);flipud(NQ(:,2))],[.8 .8 .8],'EdgeColor','none');
        hold on
    end
    plot(x,Np,'--k');
    hold on
    
    plot(x,N,'LineWidth',1.5);
    if ~tf
        hold off
    end
    
    box on
    if ischar(p.Results.covariate)
        xlabel(p.Results.covariate)
    else
        xlabel('x');
    end
    ylabel('CDF')
    legend(legendentry,'CDF(Pois)','CDF','Location','northwest')
    if strcmp(usefun,'cdf')
        ylim([0 1])
    end
end

if nargout >= 1
    varargout{1} = N(:);
    varargout{2} = x;
    varargout{3} = Np(:);
    if p.Results.confintervals || p.Results.accintervals
    varargout{4} = NQ(:,1);
    varargout{5} = NQ(:,2);
    end
    if useks
        varargout{6} = bw;
    end
end






