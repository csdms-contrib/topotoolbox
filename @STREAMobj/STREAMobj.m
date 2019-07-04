classdef STREAMobj
    
%STREAMobj Create stream object (STREAMobj)
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
%     minarea      upslope area threshold for channel initiation (default =
%                  1000)
%     unit         'pixels' (default), 'mapunits'. If you choose mapunits 
%                  provide minimum area in mapunits^2 (e.g. 1e6 m^2)
%     outlets      linear indices of drainage basin outlets (default = [])
%     channelheads linear indices of channelheads (this argument can only
%                  be used when no other parameters are set).
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
        
        narginchk(0,inf)
        
        if nargin == 0
            return
        end
        
        if ismulti(FD,true)
            error('TopoToolbox:STREAMobj','STREAMobj supports only single flow directions');
        end
        
        if nargin == 2
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
            addParamValue(p,'channelheads',[],@(x) isnumeric(x));
            
            parse(p,FD,varargin{:});
            % required
            unit    = validatestring(p.Results.unit,{'pixels', 'mapunits'});
            IX      = p.Results.outlets;
            minarea = p.Results.minarea;
            channelheads = p.Results.channelheads;
            
            
            switch unit
                case 'mapunits'
                    minarea = minarea/(FD.cellsize.^2);
            end
            
            if ~isempty(channelheads)
                W = influencemap(FD,channelheads);
            elseif ~isempty(IX) && isempty(channelheads)
                W = drainagebasins(FD,IX)>0 & flowacc(FD)>minarea;
            else
                W = flowacc(FD)>=minarea;
            end

            
            if ~any(W.Z(:))
                warning('TopoToolbox:STREAMobj',...
                    'There is no stream network that meets your criteria. \n STREAMobj returns an empty instance');
            end
            
        end
        
        % The stream obj should have only 
        % Z = false(size(W.Z));
        % Z(FD.ix) = W.Z(FD.ix);
        % Z(FD.ixc) = W.Z(FD.ixc);
        % W.Z = Z;
        
        Z = false(size(W.Z));
        Z(FD.ix)  = W.Z(FD.ix);
        I = Z(FD.ix);
        Z(FD.ixc(I)) = W.Z(FD.ixc(I));
        W.Z = Z;

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
        for r = numel(S.ix):-1:1
            distance(S.ix(r)) = distance(S.ixc(r)) + ...
                sqrt((S.x(S.ixc(r))-S.x(S.ix(r)))^2 + (S.y(S.ixc(r))-S.y(S.ix(r)))^2);
        end
    end
        
    function order = get.orderednanlist(S)
        % [dynamic property] orderednanlist returns a nan-separated vector 
        % with indices into the nodes of the STREAMobj (e.g. S.x, S.y, etc)
        
        nnal   = numel(S.x);
        nedg   = numel(S.ix);
        
        ixcix  = zeros(nnal,1);
        ixcix(S.ix) = 1:nedg;
        
        notvisited = true(nedg,1);
        row     = 1;
        counter = 1;
        
        order = nan(nedg,1);
        
        while ~isempty(row)
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
    
    function S = subgraph(S,nal)
    %SUBGRAPH extract part of the stream network
    % 
    % Syntax
    %
    %     Snew = subgraph(S,nal)
    %
    % Description
    %
    %     subgraph takes a logical node-attribute list (nal) and extracts  
    %     the nodes in the stream network where elements in nal are true.
    %
    % Input arguments
    %
    %     S       STREAMobj
    %     nal     logical node-attribute list
    %
    % Output arguments
    %
    %     Snew    STREAMobj
    %
    % Example
    %
    %     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
    %     FD = FLOWobj(DEM,'preprocess','carve');
    %     S = STREAMobj(FD,'minarea',1000);
    %     d = S.distance;
    %     nal = d>10000 & d<30000;
    %     Sn = subgraph(S,nal);
    %     plot(S)
    %     hold on
    %     plot(Sn)
    %     hold off
    %
    % See also: STREAMobj, STREAMobj/getnal, STREAMobj/modify,
    %           STREAMobj/rmnode
    % 
    % Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
    % Date: 5. June, 2019
    
    p = inputParser;
    p.FunctionName = 'STREAMobj/rmnode';
    addRequired(p,'S',@(x) isa(x,'STREAMobj'));
    addRequired(p,'nal',@(x) isnal(S,x) || isa(x,'GRIDobj'));
    parse(p,S,nal);
    
    if isa(nal,'GRIDobj')
        validatealignment(S,nal);
        nal = getnal(S,nal);
        nal = nal > 0;
    else
        nal = nal > 0;
    end
    
    if all(nal)
        % do nothing
        return
    end
    
    I = nal(S.ix) & nal(S.ixc);

    S.ix  = S.ix(I);
    S.ixc = S.ixc(I);

    IX    = cumsum(nal);

    S.ix  = IX(S.ix);
    S.ixc = IX(S.ixc);

    S.x   = S.x(nal);
    S.y   = S.y(nal);
    S.IXgrid   = S.IXgrid(nal);
    
    S     = clean(S);
     
    end
    
    function S = rmnode(S,nal)
    %RMNODE remove nodes in a stream network
    % 
    % Syntax
    %
    %     Snew = rmnode(S,nal)
    %
    % Description
    %
    %     rmnode takes a logical node-attribute list (nal) and removes the 
    %     nodes in the stream network where elements in nal are true.
    %
    % Input arguments
    %
    %     S       STREAMobj
    %     nal     logical node-attribute list
    %
    % Output arguments
    %
    %     Snew    STREAMobj
    %
    % Example
    %
    %     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
    %     FD = FLOWobj(DEM,'preprocess','carve');
    %     S = STREAMobj(FD,'minarea',1000);
    %     d = S.distance;
    %     nal = d>10000 & d<30000;
    %     Sn = rmnode(S,~nal);
    %     plot(S)
    %     hold on
    %     plot(Sn)
    %     hold off
    %
    % See also: STREAMobj, STREAMobj/getnal, STREAMobj/modify,
    %           STREAMobj/subgraph
    % 
    % Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
    % Date: 26. September, 2017    
    
    p = inputParser;
    p.FunctionName = 'STREAMobj/rmnode';
    addRequired(p,'S',@(x) isa(x,'STREAMobj'));
    addRequired(p,'nal',@(x) isnal(S,x) && islogical(x));
    parse(p,S,nal);
    S = subgraph(S,~nal);
    
    end
        
    function S = rmedge(S,eal)
    %RMEDGE remove edges in a stream network
    % 
    % Syntax
    %
    %     Snew = rmedge(S,eal)
    %
    % Description
    %
    %     RMEDGE takes a logical edge-attribute list (eal) and removes the 
    %     edges in the stream network where elements in eal are true.
    %
    % Input arguments
    %
    %     S       STREAMobj
    %     eal     logical edge-attribute list
    %
    % Output arguments
    %
    %     Snew    STREAMobj
    %
    % Example: Plot the stream network without zero gradient sections
    %
    %     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
    %     FD = FLOWobj(DEM,'preprocess','carve');
    %     S = STREAMobj(FD,'minarea',1000);    
    %     z = imposemin(S,getnal(S,DEM));
    %     I = (z(S.ix)-z(S.ixc)) == 0;
    %     Snew = rmedge(S,I);
    %     plot(Snew)
    %
    % See also: STREAMobj, STREAMobj/getnal, STREAMobj/modify
    % 
    % Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
    % Date: 26. September, 2017
    
    p = inputParser;
    p.FunctionName = 'STREAMobj/rmedge';
    addRequired(p,'S',@(x) isa(x,'STREAMobj'));
    addRequired(p,'eal',@(x) isequal(size(x),size(S.ix)) && islogical(x)); %#ok<CPROPLC>
    parse(p,S,eal);
    
    eal   = ~eal;
    
    S.ix  = S.ix(eal);
    S.ixc = S.ixc(eal);
    
    nal   = false(numel(S.IXgrid),1);
    nal(S.ix)  = true;
    nal(S.ixc) = true;

    IX    = cumsum(nal);

    S.ix  = IX(S.ix);
    S.ixc = IX(S.ixc);

    S.x   = S.x(nal);
    S.y   = S.y(nal);
    S.IXgrid   = S.IXgrid(nal);
     
    end
    
    function tf = issubgraph(S,FD)
    %ISSUBGRAPH tests if stream network is a subgraph of another stream or flow network
    %
    % Syntax
    %
    %     tf = issubgraph(S,FD)
    %     tf = issubgraph(S,S2)
    %
    % Description
    % 
    %     ISSUBGRAPH tests if a stream network S is a subgraph of the flow
    %     network in FD or another stream network in S2.
    %
    % Input arguments
    %
    %     S     STREAMobj
    %     FD    FLOWobj
    %     S2    STREAMobj
    %
    % Output arguments
    %
    %     tf    true or false scalar
    %
    % Example
    %
    %     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
    %     FD = FLOWobj(DEM,'preprocess','carve');
    %     S = STREAMobj(FD,'minarea',1000);
    %     S2 = klargestconncomps(S);
    %     % is S a subgraph of FD -> yes
    %     issubgraph(S2,FD)
    % 
    %     ans =
    % 
    %       logical
    % 
    %        1
    % 
    %     % is S a subgraph of S2 -> no
    %     issubgraph(S,S2)
    % 
    %     ans =
    % 
    %       logical
    % 
    %        0
    % 
    %     % is S2 a subgraph of S -> yes
    %     issubgraph(S2,S)
    % 
    %     ans =
    % 
    %       logical
    % 
    %        1
    % 
    % 
    % See also: STREAMobj/modify, STREAMobj/conncomps
    % 
    % Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
    % Date: 12. October, 2017
    
    
    n  = prod(S.size);
    M  = sparse(S.IXgrid(S.ix),S.IXgrid(S.ixc),true,n,n);
    if isa(FD,'STREAMobj')
        S2 = FD;
        M2 = sparse(S2.IXgrid(S2.ix),S2.IXgrid(S2.ixc),true,n,n);
        
        
    elseif isa(FD,'FLOWobj')
        M2 = FLOWobj2M(FD);
        M2 = M2>0;
    end
    tf = isequal(M2,M2|M);
    end


    function S = clean(S)

    %CLEAN remove non-connected nodes in stream networks
    %
    % Syntax
    %
    %     Sc = clean(S)
    %
    % Description
    %
    %     Modifying a STREAMobj S may sometimes generate nodes in S that have
    %     neither an incoming nor outgoing edge. This functions removes these
    %     nodes as most calculations on stream networks will dismiss them
    %     anyway.
    %
    % Input arguments
    %
    %     S     STREAMobj
    %
    % Output arguments
    %
    %     Sc    cleaned STREAMobj
    %
    % See also: STREAMobj, STREAMobj/modify
    %
    % Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
    % Date: 31. October, 2018


    % non connected nodes in the stream network are those that have neither an
    % incoming nor an outgoing edge

    M = sparse(S.ix,S.ixc,true,numel(S.x),numel(S.y));
    I = (sum(M,2) == 0) & (sum(M,1)' == 0);
    S = rmnode(S,I);
    end


end
end   
    
    
    
    