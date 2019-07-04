function [zs,exitflag,output] = crslin(S,DEM,varargin)

%CRSLIN constrained regularized smoothing of the channel length profile
%
% Syntax
%
%     zs = crslin(S,DEM)
%     zs = crslin(S,DEM,pn,pv,...)
%
% Description
%
%     Elevation values along stream networks are frequently affected by
%     large scatter, often as a result of data artifacts or errors. This
%     function returns a node attribute list of elevations calculated by
%     regularized smoothing. This function requires the Optimization 
%     Toolbox. 
%
% Input parameters
%
%     S      STREAMobj
%     DEM    Digital elevation model (GRIDobj)
%     
%     Parameter name/value pairs {default}
%
%     'K'             positive scalar that dictates the degree of stiffness
%                     {1}
%     'mingradient'   positive scalar {0}. Minimum downward gradient. 
%                     Choose carefully, because length profile may dip to
%                     steeply. Set this parameter to nan if you do not wish
%                     to have a monotonous dowstream elevation decrease.
%     'imposemin'     preprocesses the stream elevations using downstream
%                     minima imposition. {false} or true. Setting this
%                     parameter to true avoids that obviously wrong
%                     elevation values exert a strong effect on the result.
%     'nonstifftribs' {true} or false. true enables sharp concave knicks at
%                     confluences
%     'attachtomin'   true or {false}. Smoothed elevations will not exceed
%                     local minima along the downstream path. (only
%                     applicable if 'mingradient' is not nan).
%     'attachheads'   true or {false}. If true, elevations of channelheads 
%                     are fixed. (only applicable if 'mingradient' is not 
%                     nan). Note that for large K, setting attachheads to
%                     true can result in excessive smoothing and
%                     underestimation of elevation values directly 
%                     downstream to channelheads. 
%     'discardflats'  true or {false}. If true, elevations of flat sections
%                     are discarded from the analysis.
%     'knickpoints'   nx2 matrix with x and y coordinates of locations
%                     where the stiffness penalty should be relaxed (e.g. 
%                     at knickpoints, waterfalls, weirs, dams, ...).
%                     The locations should be as exact as possible but are
%                     snapped to the closest node in the stream network.
%     'precisecoords' nx3 matrix that contains the xyz coordinates of n 
%                     points whose elevations must be matched by the
%                     smoothed profile. 
%     'maxcurvature'  maximum convex curvature at any of the vertices along
%                     the profile (default: []). Setting this value to zero
%                     inhibits any convexities and the smoothed profile is 
%                     the lower concave envelope of the actual profile.
%                     This option is mainly used for the function
%                     profilesimplify, but may also be useful if the
%                     profile is strictly concave.
%     'plot'          {false} or true
%
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
%     zs  = crslin(S,DEM,'K',10,'attachtomin',true);
%     zs2 = crslin(S,DEM,'K',10);
%     plotdz(S,DEM,'color',[0.6 0.6 0.6])  
%     hold on
%     plotdz(S,zs,'color','k')
%     plotdz(S,zs2,'color','r')
%     hold off
%     legend('original data',...
%            'K=10, attached to minima',...
%            'K=10, not attached to minima')
%
% Algorithm
%
%     This algorithm uses linear nonparametric constrained regression to
%     smooth the data. The algorithm is described in Schwanghart and
%     Scherler (2017) (Eq. A6-A11).
%     
% References
%
%     Schwanghart, W., Scherler, D., 2017. Bumps in river profiles: 
%     uncertainty assessment and smoothing using quantile regression 
%     techniques. Earth Surface Dynamics, 5, 821-839. 
%     [DOI: 10.5194/esurf-5-821-2017]
%
% See also: STREAMobj/mincosthydrocon, quadprog, STREAMobj/crs, 
%           STREAMobj/quantcarve, STREAMobj/smooth
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. May, 2016


% check and parse inputs
narginchk(2,inf)

p = inputParser;
p.FunctionName = 'STREAMobj/crslin';
addParameter(p,'K',1,@(x) (isscalar(x) && x>0) || isa(x,'GRIDobj'));
addParameter(p,'imposemin',false,@(x) isscalar(x));
addParameter(p,'attachtomin',false,@(x) isscalar(x));
addParameter(p,'attachheads',false,@(x) isscalar(x));
addParameter(p,'nonstifftribs',true,@(x) isscalar(x));
addParameter(p,'mingradient',0,@(x) (isscalar(x) && (x>=0 || isnan(x))));
addParameter(p,'maxcurvature',inf);
addParameter(p,'discardflats',false);
addParameter(p,'knickpoints',[]);
addParameter(p,'precisecoords',[]);
addParameter(p,'plot',false);
addParameter(p,'weights','none',@(x) ischar(x));
addParameter(p,'nocr',false);
parse(p,varargin{:});

% get node attribute list with elevation values
if isa(DEM,'GRIDobj')
    validatealignment(S,DEM);
    z = getnal(S,DEM);
elseif isnal(S,DEM)
    z = DEM;
else
    error('Imcompatible format of second input argument')
end

if any(isnan(z))
    error('DEM or z may not contain any NaNs')
end

% get weights
wtype = lower(p.Results.weights);
switch wtype
    case 'positive'
        zimp  = imposemin(S,z);
        w = z-zimp;
        w = exp(-0.5 * w.^2);
    case 'mincost'
        zmincost = mincosthydrocon(S,z,'minmax');
        w = abs(z-zmincost);
        w = exp(-0.5 * w.^2);
    case 'interp'
        zmincost = mincosthydrocon(S,z,'interp');
        w = abs(z-zmincost);
        w = exp(-0.5 * w.^2);
    otherwise
        w = ones(size(z));
end
    
% minimum imposition
if p.Results.imposemin
    if ~exist('zimp','var')
        z = imposemin(S,z);
    else
        z = zimp;
    end
end

% upstream distance
d  = S.distance;
% nr of nodes
nr = numel(S.IXgrid);

%% Fidelity matrix
%
% This matrix is an identity matrix, since we are only interested in 
% predicting values at the node locations of the network.
if ~p.Results.discardflats
    Afid = speye(nr,nr);
    tf   = false(nr,1);
else

    % include code to allow for linear interpolation of elevations with
    % apparently wrong z values (e.g. positive residuals following carving)
    [tf,s] = identifyflats(S,imposemin(S,z));
    tf(s>0) = true;
    poi = ... %streampoi(S,'confluences','logical') | ...
          streampoi(S,'outlets','logical') | ...
          streampoi(S,'channelheads','logical');
    tf(poi) = false;
    
    % create sparse adjacency matrices weighted by the inverse distance between
    % nodes
    A = sparse(double(S.ix(tf(S.ix))),double(S.ixc(tf(S.ix))),...
        1./(d(S.ix(tf(S.ix)))-d(S.ixc(tf(S.ix)))),nr,nr); 

    A2 = sparse(double(S.ixc(tf(S.ixc))),double(S.ix(tf(S.ixc))),...
        1./(d(S.ix(tf(S.ixc)))-d(S.ixc(tf(S.ixc)))),nr,nr); 

    % since there may be more than 1 upstream neighbor, A2 is normalized such
    % that the sum of upstream weights equals the downstream weight. 
    s = 1./(sum(spones(A2),2));
    s(~tf) = 0;
    A2 = spdiags(s,0,nr,nr)*A2;

    % add A and A2
    A = A+A2;
    % and normalize again such that rows sum to one
    s = 1./(sum(A,2));
    s(~tf) = 0;
    A = spdiags(s,0,nr,nr)*A;

    % fill main diagonal with ones
    Afid  = speye(nr)-A;
    
end
        

%% Second-derivate matrix
%
% This matrix contains the finite central differences 
[I,loc] = ismember(S.ixc,S.ix);

%         i-1           i        i+1
colix  = [S.ixc(loc(I)) S.ixc(I) S.ix(I)];

% Set user-supplied knickpoints to non-stiff
if ~isempty(p.Results.knickpoints)
    xy = p.Results.knickpoints;
    [~,~,IX] = snap2stream(S,xy(:,1),xy(:,2));
    I  = ismember(S.IXgrid,IX);
    I  = I(colix(:,2));
    colix(I,:) = [];
end


val    = [2./((d(colix(:,2))-d(colix(:,1))).*(d(colix(:,3))-d(colix(:,1)))) ...
         -2./((d(colix(:,3))-d(colix(:,2))).*(d(colix(:,2))-d(colix(:,1)))) ...
          2./((d(colix(:,3))-d(colix(:,2))).*(d(colix(:,3))-d(colix(:,1))))];
      
% matrix for maximum curvature constraint      
if ~isinf(p.Results.maxcurvature)
    nrrows = size(colix,1);
    rowix  = repmat((1:nrrows)',1,3);
    Asdc   = sparse(rowix(:),colix(:),val(:),nrrows,nr);
end

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

%% Setup linear system
% balance stiffness and fidelity
% F      = p.Results.K * sqrt(size(Afid,1)/nrrows);
if ~p.Results.nocr
    F      = p.Results.K * S.cellsize^2 * sqrt(size(Afid,1)/nrrows); %norm(Afid,1)/norm(Asd,1);
    Asd    = F*Asd;
    % new with weighted scheme
    W      = spdiags(w,0,nr,nr);
    Afid   = W*Afid;
    zw     = W*z;
    C      = [Afid;Asd];
    % right hand side of equation
    b      = [zw.*(~tf);zeros(nrrows,1)]; % convex -> negative values
else
    C  = Afid;
    b  = z.*(~tf);
end

if isnan(p.Results.mingradient)
    % no minimum gradient imposition involves only solving the
    % overdetermined system of equations. 
    zs     = C\b;
    exitflag = [];
    output   = [];
        
else
    % with additional constraints we reformulate the problem for using
    % quadratic programming
    d = 1./(d(S.ix)-d(S.ixc));
    A = (sparse(S.ix,S.ixc,d,nr,nr)-sparse(S.ix,S.ix,d,nr,nr));
    e = zeros(nr,1);
    if p.Results.mingradient ~= 0
        e(S.ix) = -p.Results.mingradient;
    end
    
    if ~isinf(p.Results.maxcurvature)
        A = [A;-Asdc];
        e = [e;repmat(p.Results.maxcurvature,size(Asdc,1),1)];
    end

    % reformulate for quadprog
    H = 2*(C'*C);
    c = -2*C'*b;
    
    lb = [];
    if p.Results.attachtomin
        ub = double(z);
    else
        ub = double([]);
    end
    
    if p.Results.attachheads
        channelheads = streampoi(S,'channelheads','logical');
        Aeq = spdiags(+channelheads,0,nr,nr);
        eeq = zeros(nr,1);
        eeq(channelheads) = z(channelheads);
    else
        Aeq = [];
        eeq = [];
    end
    
    if ~isempty(p.Results.precisecoords)
        % Settings to force profiles to run through a prescribed set of 
        % elevations
        [~,~,IXp] = snap2stream(S,...
            p.Results.precisecoords(:,1),...
            p.Results.precisecoords(:,2));
        [~,locb] = ismember(IXp,S.IXgrid);
        precise = zeros(nr,1);
        precise(locb) = 1;
        
        eeqb = zeros(nr,1);
        eeqb(locb) = double(p.Results.precisecoords(:,3));
        inan = isnan(eeqb);
        eeqb(inan) = 0;
        precise(inan) = 0;
        Aeqb = spdiags(precise,0,nr,nr);
        
        if isempty(Aeq)
            Aeq = Aeqb;
            eeq = eeqb;
        else
            Aeq = Aeq+Aeqb;
            eeq = eeq + eeqb;
        end
    end
    
    % solve
    if verLessThan('optim','6.3')
        options = optimset('Display','off','Algorithm','interior-point-convex');
    else
        options = optimoptions('quadprog','Display','off');
    end
    
    % call to quadratic programming
    [zs,~,exitflag,output] = quadprog(H,c,A,e,Aeq,eeq,lb,ub,z,options);

end

if p.Results.plot
    tf = ishold(gca);
    plotdz(S,z,'color',[.6 .6 .6]);
    hold on
    plotdz(S,zs,'color','r');
    
    if ~isempty(p.Results.precisecoords)
        d = S.distance(locb);
        plot(d,p.Results.precisecoords(:,3),'s');
    end
    
    if ~tf
        hold off
    end
end


end



function [fs,s] = identifyflats(S,DEM,varargin)


if isa(DEM,'GRIDobj')
    validatealignment(S,DEM);
    z = getnal(S,DEM);
elseif isnal(S,DEM)
    z = DEM;
else
    error('Imcompatible format of second input argument')
end

fs = false(size(S.IXgrid));
s  = zeros(size(S.IXgrid));

for r=1:numel(S.ix)
    if z(S.ixc(r)) == z(S.ix(r))
%         fs(S.ix(r)) = true;
        fs(S.ixc(r)) = true;
    elseif fs(S.ix(r)) && z(S.ixc(r))< z(S.ix(r))
        s(S.ixc(r)) = 1;
    end
end

I = streampoi(S,'channelheads','logical');
fs(I) = false;
I = streampoi(S,'outlets','logical');
fs(I) = false;
s(I) = 0;
end

