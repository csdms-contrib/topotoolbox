function OUT = chiplot(S,DEM,A,varargin)

% CHI plot for bedrock river analysis
%
% Syntax
%
%     C = chiplot(S,DEM,A)
%     C = chiplot(S,DEM,A,pn,pv,...)
%
% Description
%
%     CHI plots are an alternative to slope-area plots for bedrock river
%     analysis. CHI plots are based on a transformation of the horizontal
%     coordinate that converts a steady-state river profile into a straight
%     line with a slope that is related to the ratio of the uplift rate to
%     the erodibility (Perron and Royden 2012).
%
% Input arguments
%
%     S     instance of STREAMobj. The stream network must consist of only 
%           one connected component (only one outlet may exist!)
%     DEM   digital elevation model (GRIDobj)
%     A     flow accumulation as calculated by flowacc (GRIDobj)
%
%     Parameter name/value pairs {default}
%
%     'a0': {1e6}
%     reference area in m^2
%
%     'mn':  {[]}, scalar      
%     mn is the ratio of m and n in the stream power equation. The value
%     ranges usually for bedrock rivers between 0.1 and 0.5. If empty, it
%     is automatically found by a least squares approach.
%     
%     'trunkstream': {[]}, STREAMobj
%     instance of STREAMobj that must be a subset of S, e.g. the main river
%     in the network S. The main trunk is highlighted in the plot and can
%     be used to fit the mn ratio (see pn/pv pair 'fitto').
%
%     'fitto': {'all'},'ts', 
%     choose which data should be used for fitting the mn ratio.
%     'all' fits mn to all streams in S
%     'ts' fits mn only to the trunkstream which must be provided with the
%     pn-pv pair trunkstream.
%
%     'plot': {true}, false
%     plot the CHIplot.
%
%     'mnplot': {false}, true
%     plot data for various values of mn [.1:.1:.9]
%
% Output arguments
%
%     C     structure array that contains
%     .mn       ratio of m and n
%     .beta     slope of the best fit line
%     .betase   standard error of beta
%     .a0       reference area
%     .ks       channel steepness index
%     .chi      CHI values 
%     .elev     elevation
%     .elevbl   elevation above baselevel
%     .distance distance from outlet
%     .pred     predicted elevation
%     .res      residual elevation
%
%
% Example
%
%
%
% See also: flowpathapp, STREAMobj, FLOWobj/flowacc, STREAMobj/trunk,
%           STREAMobj/modify  
%
%
% References:
%     
%     Perron, J. & Royden, L. (2013): An integral approach to bedrock river 
%     profile analysis. Earth Surface Processes and Landforms, 38, 570-576.
%     [DOI: 10.1002/esp.3302]
%     
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 13. May, 2013


% Parse Inputs
p = inputParser;         
p.FunctionName = 'chiplot';
addRequired(p,'S',@(x) isa(x,'STREAMobj'));
addRequired(p,'DEM', @(x) isa(x,'GRIDobj'));
addRequired(p,'A', @(x) isa(x,'GRIDobj'));

addParamValue(p,'mn',[],@(x) isscalar(x) || isempty(x));
addParamValue(p,'trunkstream',[],@(x) isa(x,'STREAMobj') || isempty(x));
addParamValue(p,'plot',true,@(x) isscalar(x));
addParamValue(p,'mnplot',false,@(x) isscalar(x));
addParamValue(p,'fitto','all');
addParamValue(p,'a0',1e6,@(x) isscalar(x) && isnumeric(x));
addParamValue(p,'betamethod','ls',@(x) ischar(validatestring(x,{'ls','lad'})));
addParamValue(p,'mnmethod','ls',@(x) ischar(validatestring(x,{'ls','lad'})));

parse(p,S,DEM,A,varargin{:});
S     = p.Results.S;
DEM   = p.Results.DEM;
A     = p.Results.A;
fitto = p.Results.fitto;
betamethod = validatestring(p.Results.betamethod,{'ls','lad'});
mnmethod   = validatestring(p.Results.mnmethod,{'ls','lad'});

% to which stream should the data be fitted? Trunkstream or all streams?
if ischar(fitto)
    fitto = validatestring(fitto,{'all','ts'});    
    if strcmpi(fitto,'ts') && isempty(p.Results.trunkstream);
        error('TopoToolbox:wronginput',...
             ['You must supply a trunkstream, if you use the parameter \n'...
              'fitto together with the option ts']);
    end
end

% nr of nodes in the entire stream network
nrc = numel(S.x);
M   = sparse(double(S.ix),double(S.ixc),true,nrc,nrc);
% find outlet
outlet = sum(M,2) == 0 & sum(M,1)'~=0;
if nnz(outlet)>1
    % there must not be more than one outlet (constraint could be removed
    % in the future).
    error('TopoToolbox:chiplot',...
        'The stream network must not have more than one outlet');
end

% reference drainage area
a0   = p.Results.a0; % m^2
% elevation values at nodes
zx   = double(DEM.Z(S.IXgrid));
% elevation at outlet
zb   = double(DEM.Z(S.IXgrid(outlet)));

% a is the term inside the brackets of equation 6b 
a    = double(a0./(A.Z(S.IXgrid)*(A.cellsize.^2)));

% x is the cumulative horizontal distance in upstream direction
x    = S.distance;

% use trunkstream for fitting and display
switch fitto;
    case 'all'
        SFIT = S;
        Lib = true(size(x));
    case 'ts'
        SFIT = p.Results.trunkstream;
        [Lia,Lib] = ismember(SFIT.IXgrid,S.IXgrid);
        if any(~Lia);
            error('TopoToolbox:chiplot',...
                ['The main trunk stream must be a subset of the stream network.\n'...
                 'Map a trunk stream with flowpathtool and use the STREAMobj as \n'...
                 '3rd input argument ( flowpathtool(FD,DEM,S) ).']) ;
        end
end

% find values of the ratio of m and n that generate a linear Chi plot
% uses fminsearch
if isempty( p.Results.mn );
    mn0  = 0.5; % initial value
    mn   = fminsearch(@mnfit,mn0);
else
    % or use predefined mn ratio.
    mn   = p.Results.mn;
end

% plot different values of mn
if p.Results.mnplot
    mntest = .1:.1:.9;
    cvec = jet(numel(mntest));
    figure('DefaultAxesColorOrder',cvec);
    chitest = zeros(numel(Lib),numel(mntest));
    
    for r = 1:numel(mntest);
        chitest(:,r) = netcumtrapz(x(Lib),a(Lib).^mntest(r),SFIT.ix,SFIT.ixc);
    end
    plot(chitest,zx(Lib)-zb,'x');
    xlabel('\chi [m]')
    ylabel('elevation [m]');
    title('\chi plots for different values of mn')
    legnames = cellfun(@(x) num2str(x),num2cell(mntest),'uniformoutput',false);
    legend(legnames);
end

% calculate chi
chi = netcumtrapz(x,a.^mn,S.ix,S.ixc); %*ab.^mn
% now use chi to fit beta
switch betamethod
    case 'ls'
        % least squares
        beta = chi(Lib)\(zx(Lib)-zb);
    case 'lad'
        % least absolute deviations
        beta = fminsearch(@(b) sum(abs(b*chi(Lib) - (zx(Lib)-zb))),0.0334);
end

n    = nnz(Lib);
SSE  = sum((chi(Lib)*beta - (zx(Lib)-zb)).^2);
SSZ  = sum((zx(Lib)-mean(zx(Lib))).^2);
R2   = 1-(SSE/SSZ);

betase = sqrt((SSE/(n-2))./(sum((chi(Lib)-mean(chi(Lib))).^2)));


if p.Results.plot;
    figure
    % plot results
    order = S.orderednanlist;
    I     = ~isnan(order);
    c     = nan(size(order));
    c(I)  = chi(order(I));
    zz    = nan(size(order));
    zz(I) = zx(order(I))-zb;
    
    plot(c,zz,'-','color',[.5 .5 .5]);
    hold on

    if ~isempty( p.Results.trunkstream )
        switch p.Results.fitto
            case 'all'
                % check trunkstream
        
                ST    = p.Results.trunkstream;
                [Lia,Lib] = ismember(ST.IXgrid,S.IXgrid);
        
                if any(~Lia);
                    error('TopoToolbox:chiplot',...
                        ['The main trunk stream must be a subset of the stream network.\n'...
                        'Map a trunk stream with flowpathtool and use the STREAMobj as \n'...
                        '3rd input argument ( flowpathtool(FD,DEM,S) ).']) ;
                end
            case 'ts'
                 ST = SFIT;
        end
        
        order = ST.orderednanlist; 
        I     = ~isnan(order);
        c     = nan(size(order));
        chifit = chi(Lib);
        zxfit  = zx(Lib);
        c(I)  = chifit(order(I));
        zz    = nan(size(order));
        zz(I) = zxfit(order(I))-zb;

        plot(c,zz,'k-','LineWidth',2);
    end
    
    refline(beta,0);
    hold off
    xlabel('\chi [m]')
    ylabel('elevation [m]');
end

% write to output array
if nargout == 1;
    
    OUT.mn   = mn;
    OUT.beta = beta;
    OUT.betase = betase;
    OUT.a0   = a0;
    OUT.ks   = beta*a0^mn;
    OUT.R2   = R2;
    
    [OUT.x,...
     OUT.y,...
     OUT.chi,...
     OUT.elev,...
     OUT.elevbl,...
     OUT.distance,...
     OUT.pred,...
     OUT.area] = STREAMobj2XY(S,chi,DEM,zx-zb,S.distance,beta*chi,A.*(A.cellsize^2));
     OUT.res   = OUT.elevbl-OUT.pred;

end
    


%% fitting function
function sqres = mnfit(mn)

% calculate chi with a given mn ratio
% and integrate in upstream direction
CHI = netcumtrapz(x(Lib),a(Lib).^mn,SFIT.ix,SFIT.ixc);%*ab.^mn
% normalize both variables
CHI = CHI ./ max(CHI);
z   = zx(Lib)-zb;
z   = z./max(z);
% calculate the residuals and minimize their squared sums
switch mnmethod
    case 'ls'
        sqres = sum((CHI - z).^2);
    case 'lad'
        sqres = sum(sqrt(abs(CHI-z)));
end

end

end


function z = netcumtrapz(x,y,ix,ixc)
% cumtrapz along upward direction in a directed tree network

z = zeros(size(x));
for lp = numel(ix):-1:1;
    z(ix(lp)) = z(ixc(lp)) + (y(ixc(lp))+(y(ix(lp))-y(ixc(lp)))/2) *(abs(x(ixc(lp))-x(ix(lp))));
end
end




