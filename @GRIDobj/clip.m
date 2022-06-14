function [DEM,C] = clip(DEM,C)

%CLIP clip a GRIDobj with a polygon or another GRIDobj
%
% Syntax
%
%     DEMc = clip(DEM,MS)
%     DEMc = clip(DEM,C)
%     [DEMc,Cc] = clip(...)
% 
% Description
%
%     This function clips a GRIDobj based on a polygon feature MS or
%     another GRIDobj C. C should contain a logical array. If the
%     underlying class is numeric, DEM is clipped to the non-nan pixels in
%     C. If the underlying class is integer, then DEM will be clipped to
%     the nonzero pixels in C. Pixels in DEM will be nan if they are not
%     within the clip features. If DEM contains a logical or integer
%     matrix, then pixels outside the clip features will be zero.
%
% Input arguments
%
%     DEM     GRIDobj
%     MS      mapping structure array (see polygon2GRIDobj)
%     C       GRIDobj. If C does not spatially align with DEM, then C will
%             be resampled (see GRIDobj/resample) using nearest neighbor
%             interpolation
% 
% Output arguments
%
%     DEMc    clipped GRIDobj
%     C       GRIDobj that can be used as clip feature 
%
% See also: GRIDobj/crop, GRIDobj/polygon2GRIDobj
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 27. Marchs, 2019


% Process input
if isstruct(C)
    
    C = polygon2GRIDobj(DEM,C);
elseif isa(C,'GRIDobj')
    
    % Check the data format and apply different ways to get a logical array
    if isfloat(C.Z)
        C.Z = isnan(C.Z);
    elseif isinteger(C.Z)
        C.Z = C.Z ~= 0;
    end
        
    tf = validatealignment(DEM,C);
    if ~tf
          C = resample(C,DEM,'nearest');
    end
end

% Process output
if isfloat(DEM.Z)
    DEM.Z(~C.Z) = nan;
elseif isinteger(DEM.Z)
    DEM.Z(~C.Z) = 0;
elseif islogical(DEM.Z)
    DEM.Z(~C.Z) = false;
end
        
    