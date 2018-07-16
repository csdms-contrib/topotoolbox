function [zs,kp] = knickpointfinder(S,DEM,varargin)

%KNICKPOINTFINDER find knickpoints in river profiles
%
% Syntax
%
%     [zk,kp] = knickpointfinder(S,DEM)
%     [zk,kp] = knickpointfinder(S,z)
%     [zk,kp] = knickpointfinder(...,pn,pv,...)
%
% Description
%
%     Rivers that adjust to changing base levels or have diverse
%     lithologies often feature knickzones, i.e. pronounced convex sections
%     that separate the otherwise concave equilibrium profile. The profile
%     should be monotoneously decreasing (see imposemin, quantcarve, crs).
%
%     This function extracts knickpoints, i.e. sharp convex sections in the
%     river profile. This is accomplished by an algorithm that adjusts a
%     strictly concave upward profile to an actual profile in the DEM or
%     node-attribute list z. The algorithm iteratively relaxes the
%     concavity constraint at those nodes in the river profile that have
%     the largest elevation offsets between the strictly concave and actual
%     profile until the offset falls below a user-defined tolerance.
%
%     The function returns the idealized profile zk and outputs the
%     locations of the knickpoints in a structure array kp.
%
% Input parameters
%
%     S      STREAMobj
%     DEM    Digital elevation model (GRIDobj)
%     z      node-attribute list with elevations
%     
%     Parameter name/value pairs {default}
%
%     'split'         {true} or false. If set to true, quantcarve will
%                     split the network into individual drainage basins and 
%                     process them in parallel.
%     'tol'           tolerance (scalar or node attribute list). Setting 
%                     tol to inf will return a lower concave envelope of 
%                     river profile elevations
%     'verbose'       {true} or false. Toggle verbose output to the command
%                     window
%     'plot'          {false} or true. Plots the process. Only possible if
%                     'split' is set to false.
%
% Output parameters
%
%     zk       node attribute list of river profile elevations fitted to 
%              the actual profile
%     kp       structure array 
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S);
%     zs = quantcarve(S,DEM,.5,'split',false);
%     figure
%     [zk,kp] = knickpointfinder(S,zs,'split',false,'plot',true,'tol',20);
%     hold on
%     scatter(kp.distance,kp.z,kp.dz,'sk','MarkerFaceColor','r')
%     hold off
%
% See also: STREAMobj/quantcarve, STREAMobj/crs
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 25. February, 2018

% check and parse inputs
narginchk(2,inf)

p = inputParser;
p.FunctionName = 'STREAMobj/knickpointfinder';
addParameter(p,'split',1);
addParameter(p,'knickpoints',[]);
addParameter(p,'plot',true);
addParameter(p,'tol',100);
addParameter(p,'verbose',true);
addParameter(p,'toltype','range');
parse(p,varargin{:});

if p.Results.split 
    plt = false;
    verbose = false;
else
    plt = p.Results.plot;
    verbose = p.Results.verbose;
end

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
z = double(z);

% ensure downstream decreasing profiles
z = imposemin(S,z);


%% Run in parallel ------------------------------------------------------
%  The code can run in parallel by distributing individual catchments to
%  workers. This involves some work.

if p.Results.split
    params = p.Results;
    params.split = false;
    params.plot  = false;
    [CS,locS] = STREAMobj2cell(S);
    
    tol = params.tol;
    TOL = cell(numel(CS),1);
    if ~isscalar(tol)        
        for r = 1:numel(CS)
            TOL{r} = tol(locS{r},:);
        end
    else
        [TOL{:}] = deal(tol);
    end
    params = rmfield(params,'tol');
    
    Cz  = cellfun(@(ix) z(ix),locS,'UniformOutput',false);
    Czs = cell(size(CS));
    Ckp = cell(size(CS));
    n   = numel(CS);

    parfor r = 1:n
        [Czs{r},Ckp{r}] = knickpointfinder(CS{r},Cz{r},'tol',TOL{r},params);
    end
    
    zs = nan(size(z));
    kp.n   = 0;
    kp.x   = [];
    kp.y   = [];
    kp.distance = [];
    kp.z   = [];
    kp.IXgrid  = [];
    kp.order = [];
    kp.dz = [];   
    kp.nal = false(size(z));
    
    for r = 1:numel(CS)
        zs(locS{r}) = Czs{r};
        if Ckp{r}.n > 0
        kp.nal(locS{r}) = Ckp{r}.nal;
        
        kp.n   = Ckp{r}.n+kp.n;
        kp.x   = [kp.x; Ckp{r}.x];
        kp.y   = [kp.y; Ckp{r}.y];
        kp.distance = [kp.distance; Ckp{r}.distance];
        kp.z   = [kp.z; Ckp{r}.z];
        kp.IXgrid  = [kp.IXgrid; Ckp{r}.IXgrid];
        kp.order = [kp.order; Ckp{r}.order];
        kp.dz = [kp.dz; Ckp{r}.dz];
        end

    end
    return
end

% Parallel computing goes until here. ------------------------------------

%% Knickpoint identification starts here

% upstream distance
d  = S.distance;
% nr of nodes
nr  = numel(S.IXgrid);

%% Structure array with output
if isempty(p.Results.knickpoints)
    kp.n   = 0;
    kp.x   = [];
    kp.y   = [];
    kp.distance = [];
    kp.z   = [];
    kp.IXgrid  = [];
    kp.order = [];
    kp.dz = [];
    kp.nal = [];
else
    kp    = p.Results.knickpoints;
end

%% Find knickpoints
CC = sparse(S.ix,S.ixc,true,nr,nr);
CC = speye(nr,nr) | CC | CC';

keepgoing = true;
counter   = 0;

if plt
    hh = plotdz(S,z,'color',[.6 .6 .6],'LineWidth',1.5);drawnow
    set(gcf,'color','w')
    hold on
    drawnow
    %if counter == 0
    %gif('knickpointfinder.gif','DelayTime',0.5,'LoopCount',inf,'frame',gcf)
    %end
    
end

if verbose
    disp([datestr(now) ' -- Starting']);
end

toltype = validatestring(p.Results.toltype,{'range','lowerbound'});


while keepgoing
    
    counter = counter+1;
    
    zs = lowerenv(S,z,kp.nal);

    if ~isscalar(p.Results.tol) && strcmp(toltype,'lowerbound')
        II = zs < p.Results.tol(:,1);
    else
        II = true(size(zs));
    end
    
    if counter == 1
        
        [dz,ix] = max((z-zs).*II);
        if ~isscalar(p.Results.tol)
            I = dz>0;
        else
            I = dz >= p.Results.tol;
        end
        
    else
        % calculate connected components separated by knickpoints
        CC(ix,:) = 0;
        CC(:,ix) = 0;
        CC = CC | speye(nr,nr);
        [~,pp,~,r] = dmperm(CC);
        nc = length(r) - 1;
        
        % find maximum in each connected component
        dztemp = zeros(nc,1);
        ixtemp = zeros(nc,1);
        for tt = 1:nc
            ixx = pp(r(tt):r(tt+1)-1);
            
            [dztemp(tt),ixtt] = max((z(ixx)-zs(ixx)).*II(ixx));
            ixtemp(tt) = ixx(ixtt);
        end
        
        switch p.Results.toltype
            case 'lowerbound'
                I = dztemp > 0;
            otherwise 
                if isscalar(p.Results.tol)
                    I = dztemp >= p.Results.tol;

                else
                    I = dztemp >= p.Results.tol(ixtemp);
                end
        end
        dz = dztemp(I);
        ix = ixtemp(I);       
    end
    
    if ~(any(I))
        break
        
    end
        
    IX     = S.IXgrid(ix);
    
    kp.n   = kp.n + numel(IX);
    kp.x   = [kp.x;S.x(ix)];
    kp.y   = [kp.y;S.y(ix)];
    kp.z   = [kp.z; z(ix)];
    kp.distance  = [kp.distance; d(ix)];
    kp.IXgrid    = [kp.IXgrid;IX];
    kp.order     = [kp.order; repmat(counter,numel(IX),1)];
    kp.dz        = [kp.dz;dz];
    kp.nal       = ismember(S.IXgrid,kp.IXgrid);
    
    % plot
    if plt
        if exist('hh2','var')
            delete(hh2)
        end
        hh2 = plotdz(S,zs,'color','k');
        drawnow
        %gif
    end
    if verbose
        disp([datestr(now) ' -- Iteration: ' num2str(counter) ', ' num2str(kp.n) ' knickpoints, max dz: ' num2str(max(dz))]);
    end
end

if plt
    hold off
end

end  

%% Subfunctions    
    
function z = lowerenv(S,z,ix)

% lower envelope of a channel length profile
%
% Syntax
%
%     zl = lowerenv(S,z)
%     zl = lowerenv(S,z,kn)
%
% Description
%
%     lowerenv returns the lower envelope, i.e. the lower convex hull of a
%     length profile given by the stream network S and elevation z.
%
% Input arguments
%
%     S      STREAMobj
%     z      elevation (node attribute list)
%     kn     logical vector (node attribute list) with knickpoints
%
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 10. March, 2017



kn = false(size(S.x));
if nargin == 3;
    if isempty(ix);
    elseif islogical(ix);
        kn = ix;
    else
        kn(ix) = true;
    end
end

% nals
d = distance(S);
trib = streampoi(S,'confl','logical');

nrc = numel(S.x);
ix  = S.ix;
ixc = S.ixc;

ixcix  = zeros(nrc,1);
ixcix(ix) = 1:numel(ix);

onenvelope = true(nrc,1);

for r = numel(S.ixc):-1:1;
    s  = ix(r);
    ss = ixc(r);
    
    if onenvelope(s) || trib(ss)
    IX = allpred(ix,ixc,s,nrc,kn,r);
    s  = ss;
    if isempty(IX)
        continue
    end
    gg = z(IX)-z(s);
    dd = d(IX)-d(s);
    gg = gg./dd;
    [~,ii] = sortrows([gg dd],[1 -2]); 
    IX = IX(ii(1));
    g  = gg(ii(1));
    
    t = IX;
    ixcix(ss) = 0;
    while ixcix(t) ~= 0
        t2 = ixc(ixcix(t));
        z(t2) = z(t)-g.*(d(t)-d(t2));
        onenvelope(t2) = false;
        ixcix(t) = 0;
        t  = t2;
        
    end
    end
    
end

end


    function IX = allpred(ix,ixc,s,nr,kn,startix) 
        I = false(nr,1);
        I(s) = true;
        for r = startix:-1:1
            I(ix(r)) = (I(ix(r)) | I(ixc(r))) & ~kn(ixc(r)) ;
        end
        IX = find(I);
    end

