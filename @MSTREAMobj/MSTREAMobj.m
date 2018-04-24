classdef MSTREAMobj
    
%MSTREAMobj Create multiple flow stream object (STREAMobj)
%
% Syntax
%
%     S = MSTREAMobj(FDm,INI)
%
% Description
%
%     An instance of multiple stream object (MSTREAMobj) encapsulates the
%     information on geometry and connectivity of a stream network based on
%     the multiple flow direction of a digital elevation model and a
%     logical raster that indicates the points of stream initiation.
%     MSTREAMobj provides access to various methods that investigate
%     properties of a stream network and associated data.
%
% Input arguments
%
%     FDm    instance of multiple flow direction object (FLOWobj) derived
%            by FLOWobj(DEM,'multi' or 'Dinf')
%     INI    logical grid (GRIDobj) that indicates stream locations (e.g.
%            obtained from thresholding the flow accumulation raster)
%
% Output arguments
%
%     MS     instance of MSTREAMobj
%
% Examples
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     DEM = fillsinks(DEM);
%     FDm = FLOWobj(DEM,'dinf');
%     % Randomly distribute stream initation points on all
%     % regional maxima
%     INI = GRIDobj(DEM,'logical');
%     INI.Z(545470) = true;
%     MS = MSTREAMobj(FDm,INI);
%     CS = getstream(MS,20,'method','mult');
%     imageschs(DEM,[],'colormap',[1 1 1]);
%     hold on; 
%     cellfun(@(x) plot(x,'r'),CS);
%
%  
% See also: STREAMobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 23. February, 2013
    
    
    
properties
    size      % size of instance of GRIDobj from which STREAMobj was derived
    ix        % [edge attribute] topologically sorted nodes (givers) 
    ixc       % [edge attribute] topologically sorted nodes (receivers)
    fraction  % [edge attribute] fraction delivered to pixel
    cellsize  % cellsize
    refmat    % 3-by-2 affine transformation matrix (see makerefmat)
    georef    % additional information on spatial referencing
    IXgrid    % [node attribute] linear index of stream nodes into instance of GRIDobj
    x         % [node attribute] x-coordinate vector 
    y         % [node attribute] y-coordinate vector
end


methods 
    function MS = MSTREAMobj(FD,W)
        
        narginchk(2,2)
        
        % Two input arguments: FD, W
        validatealignment(FD,W);
        p = inputParser;
        p.FunctionName = 'STREAMobj';
        addRequired(p,'FD',@(x) isa(x,'FLOWobj'));
        addRequired(p,'W', @(x) isa(x,'GRIDobj'));
        parse(p,FD,W);
        
        W = influencemap(FD,W);
        
        
        if ~any(W.Z(:))
            warning('TopoToolbox:STREAMobj',...
                'There is no stream network that meets your criteria. \n STREAMobj returns an empty instance');
        end
            
        % conn comps in W.Z must be larger than 2 pixel
        W.Z = bwareaopen(W.Z,2);

        % transfer properties from FLOWobj to STREAMobj
        MS.size     = FD.size;
        I           = W.Z(FD.ix);
        MS.ix       = double(FD.ix(I));
        MS.ixc      = double(FD.ixc(I));
        MS.fraction = double(FD.fraction(I));
        
        MS.cellsize = FD.cellsize;
        MS.refmat   = FD.refmat;
        MS.georef   = FD.georef;

        % recalculate stream indices
        IX         = zeros(FD.size,'uint32');
        IX(W.Z)    = 1:nnz(W.Z);
        MS.ix      = double(IX(MS.ix));
        MS.ixc     = double(IX(MS.ixc));
        
        I           = MS.ixc == 0;
        MS.ix(I)    = [];
        MS.ixc(I)   = [];
        
        clear IX
        MS.IXgrid  = find(W.Z);

        % get coordinate pairs
        [rows,cols] = ind2sub(MS.size,MS.IXgrid);
        xy =  [double(rows(:)) double(cols(:)) ones(numel(rows),1)] * MS.refmat;
        MS.x = xy(:,1);
        MS.y = xy(:,2);
          
    end
    
    
    function S = getstream(MS,n,varargin)
        if nargin == 1
            n = 1;
        else
            validateattributes(n,{'numeric'},{'>=',1},'getstream','n',2);
        end
        
        if n > 1
            p = gcp('nocreate');
            if isempty(p)            
                S = cellfun(@(x) getstream(MS,1,varargin{:}),num2cell(1:n),'uniformoutput',false);
            else
                S = cell(1,n);
                parfor r = 1:n
                    S{r} = getstream(MS,1,varargin{:}); %#ok<PFBNS>
                end   
            end
            return
        end
        
        p = inputParser;
        p.FunctionName = 'MSTREAMobj/getstream';
        addParameter(p,'method','random');
        parse(p,varargin{:});
        
        meth = validatestring(p.Results.method,{'random','multiplicative','additive'},3);
        
        switch meth
            case 'deterministic'
                f = MS.fraction;
            case 'uniform'
                f = rand(size(MS.fraction));
            case 'random'
                nrc = numel(MS.x);
                for r = 1:numel(MS.fraction)
                    
                end
                    
                f = rand(size(MS.fraction));
            case 'multiplicative'
                f = MS.fraction.*rand(size(MS.fraction));
            case 'additive'
                f = MS.fraction + rand(size(MS.fraction));
        end
        
        RR = (1:numel(MS.ix))';
        IX = MS.ix;
        S  = sparse(RR,IX,f,max(RR),max(IX));
        [~,ii] = max(S,[],1);
        I  = false(size(RR));
        I(ii) = true;

        ix  = MS.ix(I);
        ixc = MS.ixc(I);
        
        ch  = channelhead(MS);
        for r = 1:numel(ix)
            ch(ixc(r)) = ch(ix(r)) | ch(ixc(r));
        end
        
        I   = ch(ix) & ch(ixc);        
        
        ix  = ix(I);
        ixc = ixc(I);
        
        cc  = cumsum(double(ch));
        ix  = cc(ix);
        ixc = cc(ixc);
        
       
        
        S = STREAMobj;
        
        S.size = MS.size;
        S.ix   = ix;
        S.ixc  = ixc;
        S.cellsize = MS.cellsize;
        S.refmat = MS.refmat;
        S.georef = MS.georef;
        S.IXgrid = MS.IXgrid(ch);
        S.x     = MS.x(ch);
        S.y     = MS.y(ch);
        
    
    end
    
    function nal = channelhead(MS)
        nrc = numel(MS.x);
        M   = sparse(MS.ix,MS.ixc,true,nrc,nrc);
        nal = full(any(M,2) & ~any(M,1)');
    end
    
    function nal = outlet(MS)
        nrc = numel(MS.x);
        M   = sparse(MS.ix,MS.ixc,true,nrc,nrc);
        nal = full(any(M,1)' & ~any(M,2));
    end
    
    function h = plot(MS,varargin)
        
        xx = [MS.x(MS.ix) MS.x(MS.ixc) nan(size(MS.ixc))]';
        xx = xx(:);
        yy = [MS.y(MS.ix) MS.y(MS.ixc) nan(size(MS.ixc))]';
        yy = yy(:);
        
        ht = plot(xx,yy,varargin{:});
        
        if nargout == 1
            h = ht;
        end
    end
    
    
    function h = imagesc(MS,varargin)
        
        ix  = MS.IXgrid(MS.ix);
        ixc = MS.IXgrid(MS.ixc);
        fr  = MS.fraction;
        
        A   = zeros(MS.size);
        ch  = channelhead(MS);
        A(MS.IXgrid(ch)) = 1;
        
        for r = 1:numel(ix);
            A(ixc(r)) = A(ixc(r)) + A(ix(r)).*fr(r);
        end
        h = imagesc(sqrt(A));
        colormap(flowcolor)
    end
        
    function MS = prune(MS,minfrac)
        if nargin == 1
            minfrac = 0.05;
        end
        
        I = MS.fraction < minfrac;
        fracnew = MS.fraction;
        fracnew(I) = 0;
        
        sumfrac = accumarray(MS.ix,fracnew,size(MS.x),@sum,nan);
        isz = sumfrac == 0;
        if any(isz)
            error('The minimum fraction is too high. Choose a value less than 1/8.');
        end
        
        fracnew = fracnew./sumfrac(MS.ix);
        MS.fraction = fracnew;
        
        ch  = channelhead(MS);
        for r = 1:numel(MS.ix)
            ch(MS.ixc(r)) = (ch(MS.ix(r)) & fracnew(r)>0) | ch(MS.ixc(r));
        end
        
        ix  = MS.ix;
        ixc = MS.ixc;
        
        I   = ch(ix) & ch(ixc);    
        ix  = ix(I);
        ixc = ixc(I);
        
        cc  = cumsum(double(ch));
        ix  = cc(ix);
        ixc = cc(ixc);
        
        MS.ix  = ix;
        MS.ixc = ixc;
        MS.IXgrid = MS.IXgrid(ch);
        MS.x   = MS.x(ch);
        MS.y   = MS.y(ch);
        
        
        
        
        
    end
        
        
        
        
        
    

end
end
    
    
    
    
    