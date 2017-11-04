function OUT = slopearea(S,DEM,A,varargin)

%SLOPEAREA slope-area relation of a stream network
%
% Syntax
%
%     SA = slopearea(S,DEM,A)
%     SA = slopearea(S,z,a)
%     SA = slopearea(...,pn,pv,...)
% 
% Description
%
%     The relation between upslope area and stream gradient usually follows 
%     a power law. This function examines this relation using a stream 
%     network, a DEM (from which stream gradient is derived), and upslope 
%     area. The function aggregates gradient values into area bins and fits
%     a power law. The data (binned values, statistical parameters and
%     measures) are returned in the structure array SA.
%     
% Input arguments
%
%     S    stream network (class: STREAMobj)
%     DEM  digital elevation model (class: GRIDobj)
%     A    flow accumulation as derived from the function flowacc
%          (class: GRIDobj). 
%     z    node-attribute list (nal) with elevations
%     a    node-attribute list with flow accumulation
%     
%     parameter name/value pairs {default}
%
%     'areabins'       number of bins used to aggragate upslope area {100}.
%                      Note that the value supplied might differ from the
%                      final number of bins, since only the bins are
%                      returned that contain data.
%     'areabinlocs'    determines the way, how the area values of the bin
%                      centers are derived. Default is {'median'}, e.g.,
%                      each bin center is ultimately derived by calculating
%                      the median of the area values in each bin. 'mean' 
%                      uses the average of all area values. 'center'
%                      determines the location of the bin center by as the 
%                      0.5*(edge(i)+edge(i+1)) where edge contains the area
%                      values that limit each bin.
%     'gradaggfun'     determines how slope values are aggragated in each
%                      bin. Can be either {'mean'} or 'median'.
%     'fitmethod'      by default, this function uses a nonlinear fitting
%                      approach to avoid error bias introduced by taking
%                      the logarithms of gradient and area. The nonlinear
%                      regression lets you choose between least squares 
%                      {'ls'}, and least absolute deviations 'lad'. If you
%                      prefer a linear fit using least squares with 
%                      logtransformed data, set this option to 'logtrls'.
%     'theta'          by default []. Provide a positive value (e.g. 0.45)
%                      to only fit the scale parameter keeping the exponent
%                      at a fixed value.
%     'hist2'          if the option 'plot' is true, setting hist2 to true
%                      will result in a 2d density plot to visualize the
%                      distribution of the point cloud of the entire data.
%                      Default is {'false'}.
%     'histbins'       two element vector with number of bins used for the 
%                      2d histogram. Default ist [100 100], where the first
%                      element refers to the number of bins used for area
%                      and the second element refers to the number of bins
%                      used for counting gradient values.
%     'plot'           Plot it. Default is {true}.
%     'mingradient'    Set minimum gradient {0.0001} since fitting does not
%                      allow gradients <= 0 
%     'streamgradient' {'forward'} or 'robust'. See the parameter method of
%                      the function STREAMobj/gradient for further
%                      explanations.
%
% Output arguments
%
%     SA   structure array with following fields
%
%          .a       binned area values
%          .g       aggregated gradients
%          .ks      channel steepness index 
%          .theta   channel concavity (note that this parameter is (against
%                   convention reported as a negative value).
%          .hHist   surface handle to the histogram, if plotted
%          .hPoints line handle to the binned, empirical data, if plotted
%          .hLine   line handle to the fitted line, if plotted
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     % This DEM is not well suited to show an application of slope area
%     % plots. The DEM has various sinks located along streams resulting in
%     % many zero gradients. 
%     FD  = FLOWobj(DEM,'preprocess','c');
%     A   = flowacc(FD);
%     S   = STREAMobj(FD,A>1000);     
%     SA  = slopearea(S,DEM,A);
%
% Remark
%
%     In case you are unsatisfied by the line fitted to the data or in case
%     you require confidence intervals of the parameters, please use the
%     functions nlinfit or nlintool with following syntax. nlinfit will
%     allow you several more options for fitting and robust analysis. Note
%     that both functions require the statistics toolbox, however.
%
%     [beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(...
%                SA.a,SA.g,@(b,x) b(1)*x.^b(2),[SA.ks SA.theta]);
%
%     or 
%
%     nlintool(SA.a,SA.g,@(b,x) b(1)*x.^b(2),[SA.ks SA.theta],0.05,...
%                'area','gradient')
%
%     
%
%
% See also: slopearea, chiplot
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 19. June, 2013


narginchk(3,inf)

% Parse Inputs
p = inputParser;         
p.FunctionName = 'slopearea';

validstreamgradient = {'forward' 'centered' 'robust'};
validareabinlocs = {'center' 'median' 'mean'};
validgradaggfun  = {'mean','median'};
validfitmethods  = {'ls','lad','logtrls'};

addParamValue(p,'streamgradient','forward',@(x) ischar(validatestring(x,validstreamgradient)));
addParamValue(p,'drop',10,@(x) isscalar(x) && x>0);
addParamValue(p,'imposemin',true,@(x) isscalar(x));
addParamValue(p,'areabins',100,@(x) isscalar(x) || isempty(x));
addParamValue(p,'areabinlocs','median',@(x) ischar(validatestring(x,validareabinlocs)));
addParamValue(p,'gradaggfun','mean',@(x) ischar(validatestring(x,validgradaggfun)));
addParamValue(p,'fitmethod','ls',@(x) ischar(validatestring(x,validfitmethods)));
addParamValue(p,'fitlims',[0 inf]);
addParamValue(p,'theta',[],@(x) isscalar(x) && x>0);
addParamValue(p,'hist2',false,@(x) isscalar(x));
addParamValue(p,'plot',true,@(x) isscalar(x));
addParamValue(p,'mingradient',0.0001, @(x) isscalar(x));
addParamValue(p,'histbins',[100 100], @(x) ismember(numel(x),[1 2]) && all(x>0));

parse(p,varargin{:});
gradmeth    = validatestring(p.Results.streamgradient,validstreamgradient);
areabinlocs = validatestring(p.Results.areabinlocs,validareabinlocs);
gradaggfun  = validatestring(p.Results.gradaggfun,validgradaggfun);
fitmethod   = validatestring(p.Results.fitmethod,validfitmethods);

% get node attribute list with elevation values
if isa(DEM,'GRIDobj')
    validatealignment(S,DEM);
    z = getnal(S,DEM);
elseif isnal(S,DEM);
    z = DEM;
else
    error('Imcompatible format of second input argument')
end

% get node attribute list with elevation values
if isa(A,'GRIDobj')
    validatealignment(S,A);
    a = getnal(S,A);
elseif isnal(S,A);
    a = A;
else
    error('Imcompatible format of second input argument')
end
a = a*(S.cellsize.^2);
g = gradient(S,z,'unit','tangent',...
                   'method',gradmeth,...
                   'drop',p.Results.drop,...
                   'imposemin',p.Results.imposemin);

% evaluate               
mina = min(a);
maxa = max(a);

if p.Results.hist2;
    copya = a;
    copyg = g;
end

% bin area values
if ~isempty(p.Results.areabins);
    
    edges = logspace(log10(mina-0.1),log10(maxa+1),p.Results.areabins+1);
    [~,ix] = histc(a,edges);
    
    switch areabinlocs
        case 'mean'
            a = accumarray(ix,a,[p.Results.areabins 1],@mean,nan);
        case 'median'
            a = accumarray(ix,a,[p.Results.areabins 1],@median,nan);
        case 'center'
            a = edges(1:end-1) + diff(edges)/2;
            if a(end) == maxa;
                a(end) = nan;
            end
            a = a(:);
    end
    
    switch gradaggfun
        case 'mean'
            g = accumarray(ix,g,[p.Results.areabins 1],@(x) mean(x(~isnan(x))),nan);
        case 'median'
            g = accumarray(ix,g,[p.Results.areabins 1],@(x) median(x(~isnan(x))),nan);
       
    end
    
    I = ~isnan(a) & ~isnan(g);
    a = a(I);
    g = g(I);
    
    
    
end

g(g<=0) = p.Results.mingradient;
OUT.a = a;
OUT.g = g;

if p.Results.plot
    ax = gca;
end

% 2d Histogram
if p.Results.hist2 && p.Results.plot;
    edgesa = logspace(floor(log10(mina)),ceil(log10(maxa)),p.Results.histbins(1));
    [~,ixa] = histc(copya,edgesa);
    
    ming = min(copyg);
    ming = max(ming,p.Results.mingradient);    
    maxg = max(copyg);
    edgesg = logspace(floor(log10(ming)),ceil(log10(maxg)),p.Results.histbins(2));
    
    [~,ixg] = histc(copyg,edgesg);
    
    N   = accumarray([max(ixg(:),1) ixa(:)],1,p.Results.histbins,@sum,0);
    OUT.hHist = pcolor(ax,edgesa,edgesg,N);
    
    colormap(flipud(gray));
    shading flat
    colorbar
    hold on
end
    
% Plot dots
if p.Results.plot
    OUT.hPoints = plot(ax,a,g,'s');
    xlabel('area');
    ylabel('slope')
end


%% Fitting

% gradient 
g = max(g,p.Results.mingradient);

if ~isequal(p.Results.fitlims,[0 inf])
    I    = a>=p.Results.fitlims(1) & a<p.Results.fitlims(2);
    afit = a(I);
    gfit = g(I);
else
    afit = a;
    gfit = g;
end

if isempty(p.Results.theta) 
    % Both parameters may vary
    
        
    
    % find starting values using a least squares fit on log transformed data
    beta0 = [ones(numel(afit),1) log(afit(:))]\log(gfit);
    beta0(1) = exp(beta0(1));
    
    % fit power law S = k*A^(-mn)
    switch fitmethod
        case 'logtrls'
            beta = beta0;
        case 'ls'
            beta = fminsearch(@(beta) sum((gfit - beta(1)*afit.^beta(2)).^2),beta0);
        case 'lad'
            beta = fminsearch(@(beta) sum(abs(gfit - beta(1)*afit.^beta(2))),beta0);
    end
    
    OUT.ks = beta(1);
    OUT.theta = beta(2);

else
    % theta is fixed
    
    % slope residuals
    theta = -p.Results.theta;
    gres  = gfit./(afit.^theta);
    
    % find starting values using a least squares fit on log transformed data
    beta0 = mean(gres);
    
    % fit power law S = k*A^(-mn)
    switch fitmethod
        case 'logtrls'
            beta = beta0;
        case 'ls'
            beta = fminsearch(@(beta) sum((gfit - beta*afit.^theta).^2),beta0);
        case 'lad'
            beta = fminsearch(@(beta) sum(abs(gfit - beta*afit.^theta)),beta0);
    end
    
    OUT.ks = beta(1);
    OUT.theta = theta;
end

%% Plot
if p.Results.plot
    hold on
    aeval = logspace(log10(min(afit)),log10(max(afit)),10);
    geval = OUT.ks(1)*aeval.^OUT.theta;
    OUT.hLine = plot(ax,aeval,geval,'k-','LineWidth',1.5);
    
    set(ax,'Xscale','log','Yscale','log');
    hold off
end


         