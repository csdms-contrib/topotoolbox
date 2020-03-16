function [FD,C,S,i,k,p,A,dx_ik,kk,ii,dx_centered] = updateDrainDir(H1,BORDER,W,p,X,Y,varargin)
% Function to update the drainage network
%
% Syntax
%
%       [FD,C,S,i,k,p,DA,A,dx_ik] = updateDrainDir(H1,BORDER,W,p,X,Y,R)
%
% Description
%
%       Function to update the drainage network
%
% Input
%
%       H1        DEM (digital elevation model) (GRIDobj) 
%       BORDER    (GRIDobj) produced with getBORDER function
%       W         Weights
%       p         structure array with parameter definitions (see ttlemset)
%       X         gridded X distance Y         gridded Y distance R
%       Randomized matrix with dimensions of DEM.Z
%
%
% Output
%
%       FD        Flow Direction, FLOWobj 
%       C         River cells 
%       S         STREAMobj 
%       i         Giver 
%       k         Receiver 
%       p         updated structure array with parameter definitions (see ttlemset) 
%       A         Velocity field of advective stream power law 
%       dx_ik     distance between giver and receiver 
%       kk        Receiver of receiver 
%       ii        Giver of giver 
%       dx_centered distance between giver of giver and receiver
%
% Example
%
%
% See also:
%
% Authors:  Benjamin Campforts (benjamin.campforts@kuleuven.be)
%           Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
%
% Date:     8. July, 2016
%
%------------------------------ References -------------------------------%
% TTLEM: Campforts B., Schwanghart W., Govers G. (2016): TTLEM 1.0 : a
% numerical package for accurate simulationsng  of transient landscape
% evolution in MATLAB. Discussion paper in GMD.
%
%-------------------------------------------------------------------------%

if isempty(varargin)
    % flow direction
    if p.FlowBC
        FD = FLOWobj(H1+BORDER,'mex',true,'preprocess','carve');
    else
        FD = FLOWobj(H1,'mex',true,'preprocess','carve');
    end
else
    FD=varargin{1};
end

% upslope area ^m
DA  = flowacc(FD,W);

if p.AreaThresh > 0;
    C = DA > p.AreaThresh;
    S = STREAMobj(FD,C);
    i = S.IXgrid(S.ix);
    k = S.IXgrid(S.ixc);
else
    i = FD.ix;
    k = FD.ixc;
    S = [];
    C = GRIDobj(H1,'logical');
end

% add variable m values
if p.m_var > 0;
    %% Alternative: random error, scaled with log(DA).^2, adapted frrom Grimaldi 2005
    kc_1=1;kc=1;
    errR=(kc_1/kc^2)./(log(DA.Z).^2);
    signR=rand(size(DA.Z));
    signR(signR<=0.5)=-1;
    %Scaling
    errR=p.m_var*sign(signR).*errR./max(errR(:));
    m_var = p.m+errR*p.m;
    % figure scatter(DA.Z(:),abs(errR(:)))
else
    m_var = p.m;
end

% get velocity field
A = p.Kw.* DA.Z.^m_var;
% horizontal distance
dx_ik = sqrt((X(i)-X(k)).^2 + (Y(i)-Y(k)).^2);

%% In case of second order TVD scheme, calculate upsteream cells and second downstream cells
if strcmp(p.riverInc,'TVD_FVM')    
    [ii,kk] =get1up1down(i,k,A);
    ind_ii=ii;
    ind_ii(isnan(ii))=1;
    dx_ii = sqrt((X(ind_ii)-X(i)).^2 + (Y(ind_ii)-Y(i)).^2);
    dx_centered=(dx_ii+dx_ik)/2;
    dx_centered(isnan(ii))=dx_ik(isnan(ii));    
else
    kk=[];
    ii=[];
    dx_centered=[];
end