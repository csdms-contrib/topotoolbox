function [OUT] = mapswath(SW,M,varargin)
% obtain Z values along a SWATHobj from an arbitrary GRIDobj
%
% Syntax
%
%     [OUT] = mapswath(SW,M)
%     [OUT] = mapswath(SW,M,method)
%
% Description
%
%     MAPSWATH interpolates Z values from the GRIDobj M, at the data point
%     locations of the SWATHobj SW. If SW and M are of different spatial
%     reference (projection), MAPSWATH uses the mapping toolbox to convert
%     between units.
%
% Input arguments
%
%     SW      instance of SWATHobj
%     M       instance of GRIDobj
%     method  interpolation method, can be any of 'linear','nearest',
%             'spline','pchip','cubic'; by default 'linear'
%
% Output arguments
%
%     OUT    instance of SWATHobj
%
% Examples
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     SW = SWATHobj(DEM,'dx',200,'dy',200);
%     G = gradient8(DEM,'degree');
%     SWG = mapswath(SW,G); 
%
%
% Author: Dirk Scherler (scherler[at]caltech.edu)
% Date: June, 2013


narginchk(2,3)
if nargin == 2;
    method = 'linear';
else
    method = validatestring(varargin{1},...
        {'linear','nearest','spline','pchip','cubic'},'mapswath','method',3);
end

% check inputs
if ~isa(SW,'SWATHobj')
    error('First input hast to be of class SWATHobj');
elseif ~isa(M,'GRIDobj')
    error('Second input hast to be of class GRIDobj');
else
    
    if isempty(M.georef)
        warning('Input GRIDobj has no spatial reference.')
        OUT = SW;
    else
        
        if ~strcmp(SW.georef.Projection,M.georef.Projection) || ...
            ~isequal(SW.georef.RefMatrix,M.georef.RefMatrix)

            % inputs have different geospatial reference and/or extent
            % -> work with lat,lon
            [OUT] = swath2latlon(SW);
            switch M.georef.ModelType
                case 'ModelTypeGeographic'
                    % nothing needs to be done
                case 'ModelTypeProjected'
                    for i = 1 : length(OUT.xy0)
                        [OUT.xy0{i}(:,1),OUT.xy0{i}(:,2)] = projfwd(M.georef,OUT.xy0{i}(:,2),OUT.xy0{i}(:,1));
                        [OUT.xy{i}(:,1),OUT.xy{i}(:,2)] = projfwd(M.georef,OUT.xy{i}(:,2),OUT.xy{i}(:,1));
                        [OUT.X{i},OUT.Y{i}] = projfwd(M.georef,OUT.Y{i},OUT.X{i});
                    end
                otherwise
                    error('GRIDobj has unknown map projection. Conversion not possible.')
            end
        else
            % inputs have same geospatial reference
            OUT = SW;
        end
    end
    for i = 1 : length(OUT.xy0)
        %OUT.Z{i} = nan(size(OUT.X{i}));
        OUT.Z{i}(:) = nan;
        ix = find(~isnan(SW.Z{i}));
        OUT.Z{i}(ix) = interp(M,OUT.X{i}(ix),OUT.Y{i}(ix),method);
        % check results
        if any(isnan(OUT.Z{i}))
            warning('Output SWATHobj contains NaN. Input arguments may not fully overlap.')
        elseif all(isnan(OUT.Z{i}))
            warning('Output SWATHobj contains only NaN. Input arguments may not overlap.')
        end
    end
    
    OUT.name = {'SWATHobj created from SWATHobj:';SW.name};
    OUT.georef = M.georef;
    
end








