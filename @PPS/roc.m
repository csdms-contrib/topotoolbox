function varargout = roc(P,c,varargin)

%ROC receiver-operating characteristics of point pattern
%
% Syntax
%
%     [X,Y] = roc(P)
%     [X,Y] = roc(P,covariate)
%     [X,Y] = roc(P,mdl)
%     [X,Y] = roc(P,__,pn,pv,...)
%     [X,Y,AUC] = ...
%
% Description
%
%     roc computes the Receiver Operating Characteristic curve for a point 
%     pattern or a fitted point process modelf.
%
% Input arguments
%
%     P           instance of PPS
%     covariate   covariate of P (node attribute list or char)
%     mdl         instance of GeneralizedLinearModel as returned by
%                 fitloglinear
% 
%     Parameter name/value pairs
%     
%     plot        {true} or false
%     perfcurve   true or {false}. If true, roc uses the function
%                 perfcurve, which is part of the statistics and machine 
%                 learning toolbox.
%     
%     if perfcuve is true, than roc also takes parameter name value pairs
%     of the the function perfcurve. Specifically, these are
%
%     nboot       number of bootstrap samples for confidence intervals
%                 (Default = 20)
%     alpha       confidence level (100%*(1-alpha)) (Default = 0.05)
%
%
% Output arguments
%
%     X,Y         X and Y coordinates of an ROC curve
%     AUC         area-under-the-curve metric
%
% 
% See also: PPS, fitloglinear 
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. February, 2019


if nargin == 1
    c = P.S.distance;
end

if ~isa(c,'GeneralizedLinearModel')
    c = getcovariate(P,c);
    if size(c,2) > 1
        error('roc can handle only one covariate.')
    end
else
    c = c.Fitted.Response;
end

%% Parse inputs
p = inputParser;
p.FunctionName = 'PPS/roc';
% Add elevation
addParameter(p,'plot',true);
addParameter(p,'perfcurve',true);
addParameter(p,'NBoot',20);
addParameter(p,'BootType','norm');
addParameter(p,'Alpha',0.05);
addParameter(p,'nTvals',100);
% Parse
parse(p,varargin{:});

pp = points(P,'nal');

if p.Results.perfcurve
    TVals = quantile(c,linspace(0,1,p.Results.nTvals));
    % TVals = linspace(min(c),max(c),p.Results.nTvals);
    [X,Y,~,AUC] = perfcurve(pp,c,true,...
        'NBoot',p.Results.NBoot,...
        'BootType',p.Results.BootType,...
        'TVals',TVals);
else
    
    cp = sortrows([c pp],[1 2]);
    pp  = cp(:,2);
    
    % True positive rate
    Y = cumsum(pp);
    Y = Y./Y(end);
    
    % False positive rate
    X = cumsum(1-cp(:,2));
    X = X./X(end);
    
    AUC = trapz(X,Y);
    
    % If high densities of points are associated with low values of the
    % covariate, then we need to change the CDF
    if AUC < 0.5
        X = 1-X;
        Y = 1-Y;
    end
    
    [~,ix] = sort(X);
    X = X(ix);
    Y = Y(ix);
    AUC = trapz(X,Y);
end

% Plot results
if p.Results.plot
    tf = ishold;
    
    if p.Results.perfcurve && p.Results.NBoot > 0
        patch([X(:,2);flipud(X(:,3))],[Y(:,2); flipud(Y(:,3))],[.8 .8 .8],...
            'EdgeColor','none','FaceAlpha',0.5);
        hold on
    end
    
    plot(X(:,1),Y(:,1));
    hold on
    plot([0 1]',[0 1]','--');
    box on
    set(gca,'XLim',[0 1],'Ylim',[0 1]);
    xlabel('p')
    ylabel('roc(p)')
    if ~tf
        hold off
    end
    
end

if nargout > 0
    varargout = {X Y AUC};
end



