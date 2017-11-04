function [newz,varargout] = streamproj(S,DEM,FA,varargin)

%STREAMPROJ project stream elevations based on slope-area scaling
%
% Syntax
%
%     [zp] = streamproj(S,DEM,FA)
%     [zp] = streamproj(S,DEM,FA,pn,pv)
%     [zp,zpsig,d] = ...
%
% Description
%
%     The function streamproj allows users to project stream elevations
%     based on slope-area scaling derived from a reach along a STREAMobj.
%     This would typically be useful where longitudinal river profiles
%     display knickpoints that separate reaches with markedly different
%     channel steepness values, as resulting from the transient adjustment 
%     to a change in uplift rate, for example. Users can provide the
%     reach information in three ways. First, by providing a STREAMobj of
%     the reach, as is generated using the function modify. Second, by
%     providing the along-stream distance of the upstream and downstream
%     ends of the reach. Third, defining the reach interactively.
%     When calling the function with no parameter name-value pair, the user
%     will be prompted to interactively select a reach from a
%     distance-elevation plot of the STREAMobj S.
%     Once the reach is defined, streamproj uses the chi-formulation to
%     compute the channel steepness (ks) and concavity (theta) indices and 
%     uses these values in conjunction with Flint's law (S = ks*A^-theta)
%     to predict stream slopes based on upstream areas. It is also possible
%     to provide a theta value when calling the function (imposetheta).
%
% Input arguments
%
%     S     instance of STREAMobj
%     DEM   instance of GRIDobj with elevation
%     FA    instance of GRIDobj with upstream area
%
%
%     Parameter name/value pairs
%
%     'STREAMobj'   instance of STREAMobj
%     let's you provide a STREAMobj of the reach, which you want to base 
%     the projection on.
%
%     'distance'    scalar or [2x1] vector
%     let's you define a reach along the STREAMobj S on which the
%     projection is based on. A scalar is interpreted as minimum upstream 
%     distance and a two element vector is interpreted as minimum and 
%     maximum upstream distance.
%
%     'fitreach'    true/false
%     stream elevations are constructed downstream based on stream slopes
%     and the absolute elevation values therefore depend on the elevation
%     of the starting point. With 'fitreach' set to 'true' (default), 
%     subsequently to the reconstruction, the elevation values will be 
%     fitted to the original elevation data to reduce the influence of the 
%     starting elevation.
%
%     'imposetheta'    scalar
%     provide a theta value to constrain the concavity of the reconstructed 
%     reach.
%
%
% Output arguments
%
%     zp    projected stream elevations
%     zpsig uncertainties on the projected stream elevations as derived
%           from fitting the chi-data with a straight line
%     d     stream distance values
%
% Examples
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','c');
%     FA = flowacc(FD);
%     DEMc = imposemin(FD,DEM,0.0001);
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S,1);
%     S = trunk(S);
%     % We specify a reach along the stream using the distances of its 
%     % downstream and upstream ends 
%     [zp,zpsig,dp] = streamproj(S,DEMc,FA,'distance',[3.1e4 4.0e4]);
%     plotdz(S,DEMc), hold on
%     plot(dp,zp,'r-')
%     plot(dp,zp-zpsig,'r--')
%     plot(dp,zp+zpsig,'r--')
% 
%     % Alternatively, you can select interactively, in which case the
%     % result will be plotted automatically
%     zp = streamproj(S,DEMc,FA);
% 
%     % try out setting the theta value initially
%     zp = streamproj(S,DEMc,FA,'imposetheta',0.75);
%
%     
% See also: STREAMobj, STREAMobj/trunk, STREAMobj/modify, STREAMobj/chiplot
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: Jan 2016


narginchk(3,7)

% Parse Inputs
p = inputParser;         
p.FunctionName = 'streamproj';
addRequired(p,'S',@(x) isa(x,'STREAMobj'));
addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
addRequired(p,'FA',@(x) isa(x,'GRIDobj'));

addParameter(p,'STREAMobj',[],@(x) isa(x,'STREAMobj'));
addParameter(p,'distance',[],@(x) isnumeric(x) && numel(x)<=2);
addParameter(p,'fitreach',true,@(x) islogical(x));
addParameter(p,'imposetheta',[],@(x) isnumeric(x));


parse(p,S,DEM,FA,varargin{:});
S = p.Results.S;
R = p.Results.STREAMobj;

% Get elevation, distance, slope for the stream
[~,~,Sz,Sd,Sa] = STREAMobj2XY(S,DEM,S.distance,FA);
% note that rd is increasing while rz is decreasing, that is, when marching
% through these vectors, we effectively move downstream in the drainage.
ix = ~isnan(Sz);
Sz = Sz(ix);
Sd = Sd(ix);
Sa = Sa(ix).*FA.cellsize.*FA.cellsize;


if ~isempty(p.Results.distance)
    dist = p.Results.distance;
    if numel(dist)==1
        dist = [dist; inf];
    end
    mind = dist(1);
    maxd = dist(2);
    R = modify(S,'distance',[mind maxd]);
    
    doplot = 0;
    
elseif nargin==3 || (isempty(p.Results.distance) && isempty(R))
    
    figure
    plotdz(S,DEM,'color','k')
    
    ix = [];
    hpath = [];
    ixpath = [];
    ixend  = 0;
    title('set upper (green) reach location')
    hpstart = impoint('PositionConstraintFcn',@getnearest);
    setColor(hpstart,[0 1 0])
    idstart = addNewPositionCallback(hpstart,@drawpath);
    setPosition(hpstart,getPosition(hpstart))
    
    title('set lower (red) reach location')
    hpend   = impoint('PositionConstraintFcn',@getnearestend);
    setColor(hpend,  [1 0 0])
    idend = addNewPositionCallback(hpend,@drawpath);
    setPosition(hpend,getPosition(hpend))
    
    hptitle = title('move reach locations and press any key to extract reach');
    set(gcf,'WindowKeyPressFcn',@(k,l) uiresume);
    uiwait
    
    I = false(size(S.x));
    I(ixpath) = true;
    
    mind = Sd(ixend);
    maxd = Sd(ix);
    R = modify(S,'distance',[mind maxd]);
    
    delete(hpstart);
    delete(hpend);
    delete(hptitle);
    
    doplot = 1;
    
end

% Calculate channel steepness and concavity for the reach
if ~isempty(p.Results.imposetheta)
    C = chiplot(R,DEM,FA,'mn',p.Results.imposetheta,'plot',false);
else
    C = chiplot(R,DEM,FA,'plot',false);
end
mn = C.mn;
ks = C.ks;
ks_sig = C.betase*(C.a0^C.mn);

% Allocate space for new elevations
newz = nan(size(Sz));
newz_sig = newz;
% Move cursor to top of reach
[~,~,z] = STREAMobj2XY(R,DEM);
z0 = max(z);
ct = 1;
while Sz(ct)>z0
    ct = ct+1;
end
ct0 = ct;



% g = ks.*Sa.^(-mn);
% z = cumtrapz(S,g);
% plot(Sd,z)
%dz = fminsearch(costfun,'outletelev');




% Calculate new elevations
offset = 0;
newz(ct) = Sz(ct)-offset;
newz_sig(ct) = Sz(ct)-offset;
while ct<length(newz)
    % Calculate local slopes according to S-A-law
    thisA = Sa(ct);
    thisS = ks * thisA^(-mn);
    thisS_sig = (ks+ks_sig) * thisA^(-mn);
    thisd = Sd(ct);
    % Use local slopes to estimate dz
    ct = ct+1;
    nextd = Sd(ct);
    dz = (nextd-thisd)*thisS;
    dz_sig = (nextd-thisd)*thisS_sig;
    % Assign projected elevatioon
    newz(ct) = newz(ct-1)+dz;
    newz_sig(ct) = newz_sig(ct-1)+dz_sig;
end

if p.Results.fitreach
    zfit = newz(ct0:ct0+numel(R.x));
    origz = Sz(ct0:ct0+numel(R.x));
    costfun = @(dz) sum(((zfit+dz) - origz).^2);
    dz = fminsearch(costfun,0);
    newz = newz + dz;
    newz_sig = newz_sig + dz;
end

if nargout>1 % uncertainties
    varargout{1} = newz-newz_sig;
end
if nargout>2 % distance
    varargout{2} = Sd;
end

if doplot
    hold on
    plot(Sd,newz,'b-')
    plot(Sd,newz+newz-newz_sig,'b--')
    plot(Sd,newz_sig,'b--')
    xlims = get(gca,'XLim');
    ylims = get(gca,'YLim');
    tstr = sprintf('S = k_s A^{-theta}\n  k_s = %1.2e ± %1.2e\n  theta = %1.2f',ks,ks_sig,mn);
    text(max(xlims).*0.1,max(ylims).*0.85,tstr);
    hold off
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function posn = getnearest(pos)
    % GETNEAREST  snap to nearest location on stream network
    [~,ix] = min((Sd-pos(1)).^2 + (Sz-pos(2)).^2);
    posn= [Sd(ix) Sz(ix)];

end

function posn = getnearestend(pos)
    % GETNEAREST  snap to nearest location on stream network
    [~,ixend] = min((Sd-pos(1)).^2 + (Sz-pos(2)).^2);
    posn= [Sd(ixend) Sz(ixend)];

end

function drawpath(pos)
    % DRAWPATH   
    if ixend == 0
        ixpath = ix:numel(Sd);
    else
        ixpath = ix:ixend;
    end

    if ishandle(hpath)
        set(hpath,'Xdata',Sd(ixpath),'Ydata',Sz(ixpath));
    else
        hold on
        hpath = plot(Sd(ixpath),Sz(ixpath),'r','LineWidth',1.5);
        hold off
    end

end 


end







