function ttlem_out = ttlem(DEM,U,p,FD,varargin)

% TopoToolbox Landscape Evolution Model: TTLEM
%
% Syntax
%
%     ttlem_out = ttlem(DEM0,U,p)
%     ttlem_out = ttlem(DEM0,U,p,FD)
%
% Description
%
%     ttlem runs the TopoToolbox landscape evolution model. 
%     Parameters that control c and the numerical scheme as well as output
%     options are supplied to the function in the structure array p, which
%     is generated with ttlemset (see documentation of ttlemset).
%
% Input arguments
%
%     Z0      initial surface (digital elevation model) (GRIDobj)
%     U       grid same size as Z0 containing spatially variable uplift
%             rates [m/y] (GRIDobj)
%     p       structure array with parameter definitions (see ttlemset)
%     FD      initial flow directions (only to be used if 'DrainDir' set to
%             'fixed').
%
% Output arguments
%
%     ttlem_out.H1      final DEM
%     ttlem_out.S       final flow object, representing the river network 
%
%
% Example (documented in TTLEM_usersguide_1)
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     U = GRIDobj(DEM);
%     U.Z(2:end-1,2:end-1) = 1.5e-2;
%     p   = ttlemset;
%     p.TimeSpan = 1000000;
%     p.TimeStep = 10000;
%     p.diffScheme = 'imp_lin_sc';
%     p.Kw = 1e-4;
%     p   = ttlemset(p);
%     output = ttlem(DEM,U,p);
%     figure
%     imagesc(output.H1);
%
% *Note: for nonlinear hillslope diffusion, the semi-implicit Q-imp method
% of Taylor Perron is used in this software package. His algorithm can be
% downloaded free of charges at http://goo.gl/v1F0qs. Q-imp is further
% described in Perron (2011) doi:10.1029/2010JF001801. Although the method
% is implicit, some time constraints on the main timestep are still
% required as the method is not fully unconditionally stable for high
% uplift rates. High uplift rates may introduce gradients exceeding the the
% threshold slope in the evolving solution and hence introduce
% antidiffusion. To account for this, a conditional time step check is
% introduced and veryfied when model runs are initiated. This timestep is
% calcuatled with an empirical relation aand is set to a maxum value of
% 5000 m.
%
% See also: ttlemset
%
% References
%
% * TTLEM:
% Campforts B., Schwanghart W., Govers G.: TTLEM 1.0 : a numerical package
%           for accurate simulation of tansient landscape evolution in
%           MATLAB. Discussion paper in GMD.
%
% * TopoToolbox: Schwanghart, W. and Scherler, D.: Short Communication:
%           TopoToolbox 2 – MATLAB-based software for
%           topographic analysis and modeling in Earth surface sciences,
%           Earth Surf. Dyn., 2(1), 1–7,doi:10.5194/esurf-2-1-2014, 2014.
%           <https://www.researchgate.net/publication/259706134_Short_Communication_TopoToolbox_2_-_MATLAB-based_software_for_topographic_analysis_and_modeling_in_Earth_surface_sciences>
% 
% * TVD-FVM: Campforts, B. and Govers, G.:
%           Keeping the edge: A numerical method that avoids knickpoint
%           smearing when solving the stream power law, J. Geophys. Res.
%           Earth Surf., 120(7), 1189–1205, doi:10.1002/2014JF003376, 2015.
%           <https://www.researchgate.net/publication/279957445_Campforts_and_Govers-2015-JGR_ES_Keeping_the_edge_SI>
%
% Authors:  Benjamin Campforts (benjamin.campforts[at]kuleuven.be)
%           Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 15. June, 2016

%% Check time
dt       = p.TimeStep;
tspan    = p.TimeSpan;

if strcmp(p.diffScheme,'imp_nonlin_sc')
    parameter=10e-11;
    inc=max(max(U.Z(:)),1e-3);
    dt_max=DEM.cellsize*DEM.cellsize/inc*max(p.Sc(:))*parameter/(p.D*p.Kw);
    dt_max=min(dt_max,5000);
    nbSteps=ceil(tspan/dt_max);
    if (tspan/nbSteps)<dt
        p.TimeStep=tspan/nbSteps;
        dt       = p.TimeStep;
        disp(['Timestep is decreased to ' num2str(round(dt))...
            ' yr in order to keep the non linear solution stable']);
    end
end

rr = rem(tspan,dt);
if rr ~= 0
    tspan = tspan + (dt-rem(tspan,dt));
    warning('TopoToolbox:TTLEM',...
        ['Total simulation time was increased from %.1f to %.1f \n' ...
        'to be a multiple of the time step'],p.ts,tspan);
    p.ts   = tspan;
end
% number of iterations in the main loop
nriter = round(tspan/dt);

%% Check Boundary conditions
if strcmp(p.BC_Type,'Dirichlet_Matrix_Ini')
    p.BC_dir_value=DEM.Z;
end
%% Set DEM and Uplift values to double
DEM.Z = double(DEM.Z);
U.Z   = double(U.Z);

%% general values and parameters
%  cell size and area
dx  = DEM.cellsize;
dx2 = dx^2; %[m^2]
H0  = DEM;
nrc = numel(H0.Z);
siz = H0.size;
% initiate output
H1  = H0;
% Distance grid
[X,Y] = meshgrid((1:siz(2))*dx,(1:siz(1))*dx);

% Flow boundary conditions
if ~isempty(p.FlowBC)
    BORDER = getBORDER(DEM,p);
else
    BORDER=[];
end

%General boundary conditions
BC_indices=getBCIndices(p.BC_nbGhost,DEM.Z);%
BC.BC_indices=BC_indices;
BC.type=p.BC_Type;
BC.nbGhost=p.BC_nbGhost;
BC.dir_value=p.BC_dir_value;
if~isempty(p.BC_dir_DistSites)
    BC.BC_dir_Dist_Value=p.BC_dir_Dist_Value;
end


% Initial amount of water in each cell
W   = DEM;
W.Z = ones(W.size).*dx2;

% Weights
if ~isempty(p.NormPrecip)
    W = W*p.NormPrecip;
end

% variable drainage directions?
switch lower(p.DrainDir)
    case 'fixed'
        DrainDirVariable = false;
    otherwise
        DrainDirVariable = true;
end
if ~DrainDirVariable && exist('FD','var')    
    [FD,C,S,i,k,p,A,dx_ik,kk,ii,dx_centered] = updateDrainDir(H1,BORDER,W,p,X,Y,FD);
else
    %% River locations, needed for landsliding algorithm and diffusion
    [FD,C,S,i,k,p,A,dx_ik,kk,ii,dx_centered] = updateDrainDir(H1,BORDER,W,p,X,Y);    
end
rivLoc=C.Z(:);

%% Create diffusion matrix
switch lower(p.diffScheme)
    case {'imp_lin', 'imp_lin_sc'}
        % use Crank-Nicholson method for diffusion
        if p.D>0
            % Identity Matrix
            EYE      = speye(nrc);
            % calculate laplacian
            [ic,icd] = ixneighbors(DEM.Z,[],4);
            L        = sparse(ic,icd,1,nrc,nrc);
            L        = spdiags(sum(L,2),0,nrc,nrc) - L;
        end
end

%% Weight K
if ~isempty(p.K_weight)
    if strcmp(p.diffScheme,'imp_nonlin_sc')
        error('In case of Varying K_weight values, use a hillslope scheme different from imp_nonlin_sc in order to guarantee numerical stability');
    end
    p.Kw = p.Kw.*p.K_weight.Z;
end

%% visualize output
hGUI  = preparegui(H1);
set(hGUI.ax.time,'xlim',[0 tspan]);
drawnow

%% save directory
% does directory exist?
if ~isdir(p.resultsdir)
    mkdir(p.resultsdir)
end
resultsdir = fullfile(p.resultsdir,filesep);

%% Initiate simulation
% keep track of total time
t        = 0;
iter     = 0;

if p.steadyState
    ZPrev=DEM.Z+Inf;
    diffSS=zeros(nriter,1);
end

%% Run simulation
while iter < nriter && ~get(hGUI.stopbutton,'value');
    %% Track time
    t = t + dt;
    iter = iter + 1;
    
    %% Drainage network development
    if strcmp(p.DrainDir,'variable')
        [FD,C,S,i,k,p,A,dx_ik,kk,ii,dx_centered] = ...
            updateDrainDir(H1,BORDER,W,p,X,Y);
        rivLoc=C.Z(:);
    end
    if ~isempty(S);
        H1 = imposemin(S,H1); % Carving
    else
        H1 = imposemin(FD,H1);
    end
    
    %% Horizontal shortening
    if p.shortening
        numerics.cfl=.999;
        numerics.numMeth=p.shortening_meth;
        H1=lateralDisplacement(H1,p.short_x,p.short_y,dt, numerics,BC);
    end
    
    %% Hillslope adjustment to oversteepening
    
    switch lower(p.diffScheme)
        case {'imp_lin_sc', 'imp_nonlin_sc', 'only_sc'}
            rivEl=H1.Z(rivLoc);  %#ok<NASGU>
            H1 = excesstopography_v1(H1,'maxgradient',p.Sc*.95,'kernelsize',5,...
                'output','elevation','unit','tan','iterate',...
                true,'tol',1e-4,'maxiter',1e3);
        case {'imp_lin_sc', 'only_sc'} %#ok<MDUPC>
            %If the nonlinear diffusion scheme is applied, no slopes can exceed the threshold slope.
            H1.Z(rivLoc)=rivEl;
    end
    
    %% H1 before uplift and incision
    if strcmp(p.diffScheme,'imp_nonlin_sc')
        H_before=H1;
    end
    
    %% Vertical Uplift
    if ndims(U.Z)==2 %#ok<ISMAT>
        upl=U.Z;
    else
        upl=squeeze(U.Z(:,:,iter));
    end
    uplNoRiv=upl;
    uplNoRiv(rivLoc)=0;
    uplRiv=upl;
    uplRiv(~rivLoc)=0;
    H1 = H1 + p.rho_rs*dt.*uplNoRiv;
    
    %% River incision
    if strcmp(p.riverInc,'explicit_FDM')
        H1.Z  = funerosion_ex(p,H1.Z,dt,dx,A,i,k,dx_ik,uplRiv);
    elseif strcmp(p.riverInc,'TVD_FVM')
        %Non vectorized TVD scheme (comment and uncomment if required)
        %H1.Z  = funerosion_TVD(p,H1.Z,dt,A,i,k,dx_ik,kk,ii,dx_centered,uplRiv);
        %Vectorized TVD scheme
        H1.Z  = funerosion_TVD_vect(p,H1.Z,dt,A,i,k,dx_ik,kk,ii,dx_centered,uplRiv);
    elseif strcmp(p.riverInc,'implicit_FDM') && (p.n == 1)
        H1.Z  = funerosion_implin(p,H1.Z,dt,A, i,k,dx_ik,uplRiv);
    elseif strcmp(p.riverInc,'implicit_FDM') && (p.n ~= 1)
        if p.parallel
            if isempty(S);
                H1.Z = funerosion_impnlin_par(p,H1.Z,dt, A, FD,uplRiv);
            else
                H1.Z = funerosion_impnlin_par(p,H1.Z,dt, A, S,uplRiv);
            end
        else
            H1.Z = funerosion_impnlin(p,H1.Z,dt, A, i,k,dx_ik,uplRiv);
        end
    end
    rivEl=H1.Z(rivLoc);
    
    %% In case of non linear diffusion, do not alter the elevation of the
    % cells prior to diffusion and set H1 back to H1 before uplift and
    % icision. The forcing term is the combined impact of river incsion
    % and uplift rates and will be provided to the non linear diffusion
    % method.
    if strcmp(p.diffScheme,'imp_nonlin_sc')
        forcingTerm=(H1.Z-H_before.Z)/dt;
        H1=H_before;
    end
    
    %% Diffusion
    if p.D>0 && ~strcmp(p.diffScheme,'only_sc')
        switch lower(p.diffScheme)
            case {'imp_lin', 'imp_lin_sc'}
                H1=linearDiffusion(H1,EYE,p,dx2,C,nrc,L,dt);
            case 'imp_nonlin_sc'
                H1=nonlinDiff_Sc(H1,p, p.Sc,dt,forcingTerm,C);
                if any(H1.Z(:)<0)
                    error('Anttidiffusion: ...No - val allowed, decrease timestep!!');
                end
        end
    end
    
    %% Set river elevations to previously calcualted incision
    if strcmp(p.diffScheme,'imp_nonlin_sc') || ~p.DiffToRiv
        H1.Z(rivLoc) = rivEl;
    end
    
    %% Boundary conditions
    H1.Z=setBC(H1.Z,BC);
    if ~isempty(p.BC_dir_DistSites)
        H1.Z=BC_Disturbance(H1.Z,p.BC_dir_DistSites,BC);
    end
    
    %% Save
    if mod(iter,p.saveeach)==0;
        save([resultsdir p.fileprefix num2str(round(t)) '.mat'],'H1');
    end
    
    %% Plot
    if mod(iter,p.ploteach)==0;
        tspanstr = num2str(tspan);
        set(hGUI.hIm,'Cdata',H1.Z);
        set(hGUI.timedisp,'String',[num2str(round(t)) ' / ' tspanstr ' years (' num2str(t/tspan*100,'%3.1f') '%)'] );
        
        clim  = get(hGUI.ax.main,'clim');
        tt    = get(hGUI.tsmin,'Xdata');
        ydmin = get(hGUI.tsmin,'Ydata');
        ydmax = get(hGUI.tsmax,'Ydata');
        ydmean= get(hGUI.tsmean,'Ydata');
        tt    = [tt(:); t];
        set(hGUI.tsmin,'XData',tt,'YData',[ydmin(:);clim(1)]);
        set(hGUI.tsmax,'XData',tt,'YData',[ydmax(:);clim(2)]);
        set(hGUI.tsmean,'XData',tt,'YData',[ydmean(:);mean(H1.Z(~isnan(H1.Z)))]);
        drawnow
    end
    
    %% Check for steadyState
    if p.steadyState
        diffSS(iter)=abs(sum(ZPrev(:)-H1.Z(:)));
        if diffSS(iter)<p.SS_Value
            disp(['Steady state achieved after ' num2str(round(iter*p.TimeStep*1e-3)*1e-3) ' Myr'])
            % generate output
            ttlem_out.H1   = H1;
            ttlem_out.S   = S;
            ttlem_out.diffSS=diffSS;
            ttlem_out.SSTime=t;
            return
        end
        ZPrev=H1.Z;
    end
    
end

%% Generate output
ttlem_out.H1   = H1;
ttlem_out.S   = S;
if p.steadyState
    ttlem_out.diffSS=diffSS;
    ttlem_out.SSTime=t;
end
end


