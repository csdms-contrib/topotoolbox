function OUT = chopswath(SW,varargin)
% chop a SWATHobj into equally long segments
%
% Syntax
%
%    OUT = chopswath(SW,d)
%
%
% Description
%     
%     CHOPSWATH(SW,d) chops the SWATHobj SW into segments of length 'd'.
%     The newly created segments are successively obtained from the
%     previous segments of SW, so that towards the end of the previous
%     segments there may be newly created segments of shorter length than
%     'd'.
%
%
% Input arguments
%
%     SW     instance of SWATHobj
%     d      chopping distance (scalar)
%
% Output
%     
%     OUT    new instance of SWATHobj
%
% 
% Example 
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     FA  = flowacc(FD);
%     S = STREAMobj(FD,FA>1000);
%     S = klargestconncomps(S,1);
%     S = removeshortstreams(S,1e3);
%     SW = SWATHobj(DEM,S,'smooth',300,'plot',false);
%     SW = tidyswath(SW,FD,'both');
%     SW2 = chopswath(SW,1e3);
%     plot(SW2)
%
%
% Author: Dirk Scherler (scherler[at]caltech.edu)
% Date: Novemver, 2013



narginchk(2,2)


if ~isa(SW,'SWATHobj')
    error('First argument has to be of class SWATHobj');
end


chopdist = varargin{1};
IX = cell(size(SW.Z));
for i = 1 : length(SW.Z)
    d = SW.distx{i};
    ix = mod(d-min(d),chopdist);
    isix = [0; diff(ix)<0];
    ix1 = [1;find(isix)];
    ix2 = [find(isix)-1;length(isix)];
    IX{i} = [ones(size(ix1)).*i,ix1,ix2];
end

% Assemble new SWATHobj
OUT = SW;
OUT.xy0(:) = [];
OUT.zd0(:) = [];
OUT.smooth(:) = [];
OUT.xy(:) = [];
OUT.distx(:) = [];
OUT.disty(:) = [];
OUT.X(:) = [];
OUT.Y(:) = [];
OUT.Z(:) = [];

ct = 0;
for i = 1 : length(IX)
    for j = 1 : size(IX{i},1)
        ct = ct+1;
        kx = IX{i}(j,1);
        ix1 = IX{i}(j,2);
        ix2 = IX{i}(j,3);
        OUT.smooth{ct} = SW.smooth{kx};
        OUT.xy{ct} = SW.xy{kx}(ix1:ix2,:);
        OUT.distx{ct} = SW.distx{kx}(ix1:ix2);
        OUT.disty{ct} = SW.disty{kx}(:);
        OUT.X{ct} = SW.X{kx}(:,ix1:ix2);
        OUT.Y{ct} = SW.Y{kx}(:,ix1:ix2);
        OUT.Z{ct} = SW.Z{kx}(:,ix1:ix2);
        d0 = SW.zd0{IX{i}(j,1)}(:,2);
        ix = find(d0>=OUT.distx{ct}(1) & d0<=OUT.distx{ct}(end));
        if isempty(ix); ix = length(d0); end
        OUT.xy0{ct} = SW.xy0{kx}(ix,:);
        OUT.zd0{ct} = SW.zd0{kx}(ix,:);
    end
end








