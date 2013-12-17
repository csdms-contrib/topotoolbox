classdef STREAMobj
    
% Create stream object (STREAMobj)
%
% Syntax
%
%     S = STREAMobj(FD,W)
%     S = STREAMobj(FD,pn,pv)
%
% Description
%
%     An instance of stream object encapsulate the information on geometry
%     and connectivity of a stream network based on the flow direction of a
%     digital elevation model and a logical raster that indicates the
%     position of streams. STREAMobj provides access to various methods
%     that investigate properties of a stream network and associated data.
%
% Input arguments
%
%     FD     instance of flow direction object (FLOWobj)
%     W      logical grid (GRIDobj) that indicates stream locations (e.g.
%            obtained from thresholding the flow accumulation raster)
%
% Parameter name/value pairs
%
%     minarea    upslope area threshold for channel initiation (default =
%                1000)
%     unit       'pixels' (default), 'mapunits'. If you choose mapunits 
%                provide minimum area in mapunits^2 (e.g. 1e6 m^2)
%     outlets    linear indices of drainage basin outlets (default = [])
%
% Output arguments
%
%     S      instance of STREAMobj
%
% Examples
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','c');
%     S = STREAMobj(FD,flowacc(FD)>1000);
%     plot(S)
%
% Note: You can use flowpathapp to manually map channelheads and export an
%       instance of STREAMobj.
%  
% See also: STREAMobj/modify, FLOWobj, GRIDobj, flowpathapp
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 23. February, 2013
    
    
    
properties
    size      % size of instance of GRIDobj from which STREAMobj was derived
    ix        % [edge attribute] topologically sorted nodes (givers) 
    ixc       % [edge attribute] topologically sorted nodes (receivers)
    cellsize  % cellsize
    refmat    % 3-by-2 affine transformation matrix (see makerefmat)
    georef    % additional information on spatial referencing
    IXgrid    % [node attribute] linear index of stream nodes into instance of GRIDobj
    x         % [node attribute] x-coordinate vector 
    y         % [node attribute] y-coordinate vector
end

properties (Dependent = true)    
    distance  % [node attribute] distance from outlet (dynamic property)
    orderednanlist % nan-separated index into node attributes
end


methods 
    function S = STREAMobj(FD,varargin)
        
        narginchk(2,inf)
        
        if nargin == 2;
            % Two input arguments: FD, W
            W = varargin{1};
            validatealignment(FD,W);
            p = inputParser;         
            p.FunctionName = 'STREAMobj';
            addRequired(p,'FD',@(x) isa(x,'FLOWobj'));
            addRequired(p,'W', @(x) isa(x,'GRIDobj'));
            parse(p,FD,W);
        else 
            % many input arguments: FD, pn-pv pairs
            p = inputParser;
            p.FunctionName = 'STREAMobj';
            addRequired(p,'FD',@(x) isa(x,'FLOWobj'));
            addParamValue(p,'minarea',1000,@(x) isscalar(x) && x>=0);
            addParamValue(p,'unit','pixels',@(x) ischar(validatestring(x, ...
                            {'pixels', 'mapunits'}))); 
            addParamValue(p,'outlets',[],@(x) isnumeric(x));
            
            parse(p,FD,varargin{:});
            % required
            unit = validatestring(p.Results.unit,{'pixels', 'mapunits'});
            IX   = p.Results.outlets;
            minarea = p.Results.minarea;
            
            switch unit
                case 'mapunits';
                    minarea = minarea/(FD.cellsize.^2);
            end
            
            if ~isempty(IX)
                W = drainagebasins(FD,IX)>0 & flowacc(FD)>minarea;
            else
                W = flowacc(FD)>=minarea;
            end

            
            if ~any(W.Z(:))
                warning('TopoToolbox:STREAMobj',...
                    'There is no stream network that meets your criteria. \n STREAMobj returns an empty instance');
            end
            
        end
      
        % conn comps in W.Z must be larger than 2 pixel
        W.Z = bwareaopen(W.Z,2);

        % transfer properties from FLOWobj to STREAMobj
        S.size     = FD.size;
        I          = W.Z(FD.ix);
        S.ix       = double(FD.ix(I));
        S.ixc      = double(FD.ixc(I));
        
        S.cellsize = FD.cellsize;
        S.refmat   = FD.refmat;
        S.georef   = FD.georef;

        % recalculate stream indices
        IX        = zeros(FD.size,'uint32');
        IX(W.Z)   = 1:nnz(W.Z);
        S.ix      = double(IX(S.ix));
        S.ixc     = double(IX(S.ixc));
        
        I          = S.ixc == 0;
        S.ix(I)    = [];
        S.ixc(I)   = [];
        
        clear IX
        S.IXgrid  = find(W.Z);

        % get coordinate pairs
        [rows,cols] = ind2sub(S.size,S.IXgrid);
        xy =  [double(rows(:)) double(cols(:)) ones(numel(rows),1)] * S.refmat;
        S.x = xy(:,1);
        S.y = xy(:,2);
          
    end
    
    
    
    
    
    
    function distance = get.distance(S)
        % [dynamic property] distance from outlet
        distance = zeros(numel(S.x),1);
        for r = numel(S.ix):-1:1;
            distance(S.ix(r)) = distance(S.ixc(r)) + ...
                sqrt((S.x(S.ixc(r))-S.x(S.ix(r)))^2 + (S.y(S.ixc(r))-S.y(S.ix(r)))^2);
        end
    end
        
    function order = get.orderednanlist(S)
        % [dynamic property] orderednanlist returns a nan-separated vector 
        % with indices into the nodes of the STREAMobj (e.g. S.x, S.y, etc)
        ixcix  = zeros(size(S.x));
        ixcix(S.ix) = 1:numel(S.ix);
        
        notvisited = true(size(S.ix));
        row     = 1;
        counter = 1;
        
        order = nan(size(S.ix));
        
        while ~isempty(row);
            % initiate new stream
            order(counter)  = S.ix(row);
            % increase counter
            counter         = counter+1;
            
            while row~=0  
                % follow stream
                order(counter)   = S.ixc(row);
                % inidicate as visited
                ixcix(S.ix(row)) = 0;
                notvisited(row)  = false;
                % jump to next row
                row              = ixcix(S.ixc(row));
                % increase counter
                counter          = counter + 1;
            end
            % you have reached the end of the stream
            % separate streams with a nan
            order(counter) = nan;
            % increase counter
            counter        = counter+1;
            % find the next stream head
            row            = find(notvisited,1,'first');
        end
        
        order = order(:);
    end


end
end
    
    
    
    
    