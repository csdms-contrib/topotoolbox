function OUT = modify(SW,varargin)
% set fields in a SWATHobj to NaN according to different rules
%
% Syntax
%
%     OUT = modify(SW,'pn','pv',...)
%
% Description
%
%     MODIFY uses optional rules specified by the user to set fields in the
%     SWATHobj SW to NaN.
%
% Input arguments
%
%     SW      instance of SWATHobj
%
%     Parameter name/value pairs
%
%     'distance'   1x2 real vector
%     limits length of SWATHobj by imposing constraints on the minimum and 
%     maximum distance along swath. default is 
%
%     'width'   scalar
%     limits width of SWATHobj
%
%     'gap'   scalar
%     changes central gap of SWATHobj
%
%     'absolutez'  1x2 real vector
%     excludes data points depending on their absolute z value
%
%     'relativez'  1x2 real vector
%     excludes data points at equal distance along the SWATHobj, depending
%     on their relative z value with respect to the minimum z value at 
%     that distance
%
%     'logical'    GRIDobj
%     uses a logical GRIDobj to exclude data points along SWATHobj
%
% Output arguments
%
%     OUT    instance of SWATHobj
%
% Examples
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     FA = flowacc(FD);
%     S = STREAMobj(FD,FA>1e6/FA.cellsize/FA.cellsize);
%     S = trunk(klargestconncomps(S));
%     SW = SWATHobj(DEM,S,'dx',100,'dy',100,'width',3000,'smooth',11,'plot',false);
%     SW2 = modify(SW,'relativez',[0 100]);
%     imagesc(DEM), hold on
%     plot(SW2,'outline',false)
%
%
% Author: Dirk Scherler (scherler[at]caltech.edu)
% Date: June, 2013


OUT = SW;

% Parse inputs
p = inputParser;
p.FunctionName = 'modify';
addRequired(p,'SW',@(x) isa(x,'SWATHobj'));
addParamValue(p,'distance',[min(cell2mat(SW.distx')), max(cell2mat(SW.distx'))],@(x) isnumeric(x))
addParamValue(p,'width', max(cell2mat(SW.disty')),@(x) isnumeric(x))
addParamValue(p,'gap',min(abs(cell2mat(SW.disty'))),@(x) isnumeric(x))
addParamValue(p,'relativez',[-inf,inf],@(x) isnumeric(x))
addParamValue(p,'absolutez',[-inf,inf],@(x) isnumeric(x))
addParamValue(p,'logical',[],@(x) isa(x,'GRIDobj'))
parse(p,SW,varargin{:});

% parameters
newdistx     = p.Results.distance;
newwidth     = p.Results.width;
newgap       = p.Results.gap;
relz         = p.Results.relativez;
absz         = p.Results.absolutez;
M            = p.Results.logical;

if ~isempty(M)
    [tSW] = mapswath(SW,M,'nearest');
end


for i = 1 : length(SW.xy0)
    
    % distx
    IX = find(SW.distx{i}<newdistx(1) | SW.distx{i}>newdistx(2));
    OUT.distx{i}(IX) = nan;
    OUT.xy{i}(IX,:) = nan;
    OUT.X{i}(:,IX) = nan;
    OUT.Y{i}(:,IX) = nan;
    OUT.Z{i}(:,IX) = nan;
    
    % width & gap
    IX = find(abs(SW.disty{i})<newgap | abs(SW.disty{i})>newwidth);
    OUT.disty{i}(IX) = nan;
    OUT.X{i}(IX,:) = nan;
    OUT.Y{i}(IX,:) = nan;
    OUT.Z{i}(IX,:) = nan;
    
    minz = min(OUT.Z{i},[],1);
    tZ = OUT.Z{i}-repmat(minz,length(OUT.disty{i}),1);
    IX = tZ<relz(1) | tZ>relz(2);
    OUT.X{i}(IX) = nan;
    OUT.Y{i}(IX) = nan;
    OUT.Z{i}(IX) = nan;
    
    IX = OUT.Z{i}<absz(1) | OUT.Z{i}>absz(2);
    OUT.X{i}(IX) = nan;
    OUT.Y{i}(IX) = nan;
    OUT.Z{i}(IX) = nan;
    
    if ~isempty(M)
        IX = find(tSW.Z{i}==0);
        OUT.X{i}(IX) = nan;
        OUT.Y{i}(IX) = nan;
        OUT.Z{i}(IX) = nan;
    end
    
end % loop over SWATHobj entries


