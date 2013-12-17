function [OUT,varargout] = drainagebasins(FD,varargin)

% drainage basin delineation/catchments
%
% Syntax
%
%     L = drainagebasins(FD)
%     L = drainagebasins(FD,IX)
%     L = drainagebasins(FD,x,y)
%     L = drainagebasins(FD,SO,order)
%     [L,outlets] = ...
%     [L,x,y] = ...
%
% Description
%
%     drainagebasins determines drainage basin affiliation of cells in 
%     a digital elevation model (DEM) based on flow direction (FLOWobj).
%     When drainagebasins is supplied with the flow direction matrix, all
%     cells have a clear affiliation to a drainage basin.
%
%     drainagebasins(FD) calculates the drainage basins from all drainage
%     network outlets.
%
%     drainagebasins(FD,IX) and drainagebasins(FD,x,y) lets you specify
%     outlets indicated by the linear index vector (IX) or coordinate
%     vectors (x,y) to derive drainage basins.
%
%     drainagebasins(FD,SO,order) takes a streamorder grid (SO) calculated
%     by the function FLOWobj/streamorder and a scalar that indicates the
%     stream order for which drainage basins are to be derived. order = 1 
%     for instance, outputs all first order drainage basins. Note that
%     drainage basins are not calculated for stream links of specified 
%     order that are attached to the grids edges.
%
% Input
%
%     FD        flow direction object (FlowDirObj)
%     IX        linear index indicating position of user defined drainage 
%               basin outlets
%     x,y       coordinate pairs of user defined drainage basin outlets
%     S         stream order raster as produced by the function
%               FLOWobj/streamorder (GRIDobj)
%     order     stream order (scalar integer)
% 
% Output
%
%     L         drainage basin grid (GRIDobj, uint32)
%     outlets   linear indices of drainage basin outlet cells
%     x,y       coordinate pairs of drainage basin outlet cells
%
% Example 1
%
%     % get all drainagebasins
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','c');
%     D   = drainagebasins(FD);
%     imageschs(DEM,D)
%
% Example 2
%
%     % get all drainage basins with streamorder 3
%     S = streamorder(FD,flowacc(FD)>100);
%     D = drainagebasins(FD,S,3);
%     imageschs(DEM,D)
%
%
%
% See also: FLOWobj, GRIDobj/snap2stream, FLOWobj/streamorder,
%           GRIDobj/shufflelabel
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 10. February, 2012


narginchk(1,3);

if strcmpi(FD.type,'multi');
    error('TopoToolbox:drainagebasins',...
        'drainage basins are not defined for multiple flow directions');
end

if exist(['drainagebasins_mex.' mexext],'file')==3 && ...
   nargout==1 && nargin == 1;
    % Use mex-file
    D = drainagebasins_mex(FD.ix,FD.ixc,FD.size);
    
elseif nargin == 1;
    % Don't use mex-file
    DBcounter = 0;
    D = zeros(FD.size,'uint32');

    for r = numel(FD.ix):-1:1;
        if D(FD.ixc(r)) == 0;
            DBcounter = DBcounter+1;
            D(FD.ixc(r)) = DBcounter;
            outlets(DBcounter) = FD.ixc(r);
        end
        D(FD.ix(r)) = D(FD.ixc(r));
    end
elseif nargin > 1;
    % ix,x and y, or Stream order grid and stream order are supplied
    if nargin == 2;
        % IX is supplied
        IX = varargin{1}; 
        validateattributes(IX,{'numeric'},{'vector','integer'})
    else
        if isa(varargin{1},'GRIDobj')
            S = varargin{1};
            % find seed pixels for specific stream order
            S.Z(S.Z ~= varargin{2}) = 0;            
            I = (S.Z(FD.ix) & ~S.Z(FD.ixc));
            IX = FD.ix(I);
        else
            % SEED pixels are supplied as coordinate pairs
            IX   = coord2ind(FD,varargin{1},varargin{2});
        end
    end
 
    D = zeros(FD.size,'uint32');
    D(IX) = cast(1:numel(IX),'uint32');
        
    for r = numel(FD.ix):-1:1;
        if D(FD.ixc(r)) ~= 0 && D(FD.ix(r))==0;
            D(FD.ix(r)) = D(FD.ixc(r));
        end
    end
    
    outlets = IX;
end

%% Prepare Output
OUT = copy2GRIDobj(FD);
% write output to GRIDobj
OUT.Z = D;
OUT.zunit = '';
OUT.name  = 'drainage basins';

if nargout == 2;
    varargout{1} = outlets;
elseif nargout == 3;
    [x,y] = ind2coord(S,outlets);
    varargout{1} = x;
    varargout{2} = y;
end


    