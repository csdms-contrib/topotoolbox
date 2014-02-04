classdef FLOWobj
    
% Create flow direction object
%
% Syntax
%
%     FD = FLOWobj(DEM)
%     FD = FLOWobj(DEM,'pn','pv',...)
%     FD = FLOWobj(M,'cellsize',cs,'size',siz)
%
% Description
%
%     FLOWobj creates a flow direction object which can be used to
%     calculate terrain attributes such as flow accumulation, drainage
%     basin delineation, flow path extraction, etc. In contrast to the M
%     matrix returned by the function flowdir, the flow direction object
%     requires much less memory and the same routines are usually evaluated
%     much faster.
%
%     FD can be calculated from the original digital elevation model or
%     converted from an existing flow direction matrix. FD contains a
%     topological ordering of the DEM nodes.
%
% Input arguments
%
%     DEM    digital elevation model (Class: GRIDobj)
%     M      flow direction matrix as returned by flowdir, flowdir_single
%            (function of TopoToolbox 1)
%     
% Parameter name/value pairs   {default}
%
%  Applicable only, if calculated from GRIDobj
%
%     'preprocess'   {'fill'}, 'carve', 'none'
%            set DEM preprocessing that determines flow behavior in
%            topographic depressions and flat areas 
%     'sinks'  logical matrix same size as dem
%            true values in sinks are treated as sinks in the digital
%            elevation model and are not filled or carved, if the
%            preprocessing option fill or carve are chosen.
%     'internaldrainage' {false} or true
%            set this parameter value to true if flow directions should be
%            derived in the lowest, flat regions of internal drainage 
%            basins. By default, this parameter is set to false since this
%            information is usually not required and flow paths will stop
%            when entering flat, internally drained sections.
%     'cweight'   scalar {1}
%            adjust cost matrix if preprocessing option 'carve' has been
%            chosen. 
%     'verbose'   {false},true
%            verbose output in the command window to track computational
%            progress. Particularly interesting when working with very
%            large matrices.
%     'mex'   {'false'},true
%            controls if the mex routines should be called. If true, the mex
%            functions (see function compilemexfiles) must be available on
%            the search path. Using mex functions can increase the speed at
%            which an instance of FLOWobj is constructed.
%
% Applicable only, if calculated from flow direction matrix M
%
%     'cellsize' {1}
%            cellsize in x and y direction (scalar double) 
%     'size'   (1x2 array, e.g. size(dem))
%            the size of the DEM of which the flow direction matrix was 
%            derived from. Only applicable when function called with flow 
%            direction matrix.
%     'refmat' referencing matrix as derived from GRIDobj property refmat
%
% Output
%
%     FD     flow direction object (FLOWobj)
%
% Example
%
%     % create and evaluate flow direction object for a DEM
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     A  = flowacc(FD);
%     imageschs(DEM,log(A));
%
%
% 
% See also: GRIDobj, STREAMobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 5. January, 2013

properties(GetAccess = 'public', SetAccess = 'private')
    size      % size of instance of GRIDobj from which STREAMobj was derived
    type      % flow direction type (single, multi)
    ix        % [edge attribute] topologically sorted nodes (givers)
    ixc       % [edge attribute] topologically sorted nodes (receivers)
    fraction  % [edge attribute] fraction transfered between nodes
    cellsize  % cellsize of the grid (scalar)
    refmat    % 3-by-2 affine transformation matrix (see makerefmat)
    
end

properties(GetAccess = 'public', SetAccess = 'public')
    fastindexing = false; % set to true to initiate fast indexing
    ixcix     % indexing matrix for fast indexing
    georef    % additional information on spatial referencing
end

methods
    function FD = FLOWobj(DEM,varargin)

        if nargin == 0;
            
        else

            % Parse inputs
            p = inputParser;
            p.FunctionName = 'FLOWobj';
            expectedPreProcess = {'none','fill','carve'};

            addRequired(p,'DEM',@(x) issparse(x) || isa(x,'GRIDobj'));
            
            addParamValue(p,'size',[],@(x) isempty(x) || numel(x)==2);
            addParamValue(p,'cellsize',1,@(x) isscalar(x));
            addParamValue(p,'refmat',[]);
            addParamValue(p,'preprocess','none',@(x) any(validatestring(x,expectedPreProcess)));
            addParamValue(p,'tweight',2,@(x) isscalar(x));
            addParamValue(p,'cweight',1,@isnumeric);
            addParamValue(p,'sinks',[],@(x) isa(x,'GRIDobj'));
            addParamValue(p,'verbose',false,@(x) isscalar(x) && islogical(x));
            addParamValue(p,'streams',1);
            addParamValue(p,'weights',[]);
            addParamValue(p,'internaldrainage',false,@(x) isscalar(x) && islogical(x));
            addParamValue(p,'mex',false,@(x) isscalar(x) && islogical(x));

            parse(p,DEM,varargin{:});

            % required
            DEM        = p.Results.DEM;
            siz        = p.Results.size;
            tweight    = p.Results.tweight;
            preprocess = validatestring(p.Results.preprocess,expectedPreProcess);
            cweight    = p.Results.cweight;
            sinks      = p.Results.sinks;
            verbose    = p.Results.verbose;
            streams    = p.Results.streams;
            weights    = p.Results.weights;
            mexx       = p.Results.mex;
            
            % % % %
            if mexx
                % check if compiled mex files are really available
                ext = mexext;
                mexx = exist(['steepestneighbor_mex.' ext],'file')==3 && ...
                       exist(['tsort_mex.' ext],'file')==3;
                if mexx==0;
                    warning('TopoToolbox:FLOWobj',...
                            ['cannot find compiled mex functions on the search path. \n'...
                             'FLOWobj continues with using slower m-functions']);
                end
            end
                

            % % % % % %
            if issparse(DEM)
                % construct flow direction object from M
                M = DEM;
                FD.cellsize = p.Results.cellsize;
                FD.refmat   = p.Results.refmat;
                
                if ~isempty(siz);
                    FD.size = siz;
                else
                    fprintf(['\n       The size of the DEM was not provided. It will be \n'...
                        '       determined automatically. For large datasets it \n'...
                        '       may be more efficient to provide the size as \n'...
                        '       optional input argument.\n \n']);


                    [ic,icd] = find(M);
                    offsets  = unique(ic-icd);
                    FD.size(1) = max(abs(offsets)-1);
                    FD.size(2) = size(M,1)./FD.size(1); %#ok<CPROP>

                    fprintf('       The size of the DEM was determined to be [%d %d].\n \n',...
                        FD.size(1),FD.size(2));
                end

                mtf = ismulti(M);

                if mtf
                    [FD.ix,FD.ixc,FD.fraction] = tsort(M);
                    FD.type = 'multi';
                else
                    [FD.ix,FD.ixc] = tsort(M);
                    FD.type = 'single';
                end
                
                FD.ix = uint32(FD.ix);
                FD.ixc = uint32(FD.ixc);
                
                

            %%  construct flow direction object from DEM
            else

                FD.cellsize = DEM.cellsize;
                FD.refmat   = DEM.refmat;
                FD.georef   = DEM.georef;
                FD.size     = DEM.size;
                FD.type     = 'single';
                nrc         = numel(DEM.Z);

                switch lower(preprocess)                   
                    
                    case 'fill'
                        if verbose
                            disp([datestr(clock) ' -- Sink filling'])
                        end
                        
                        if isempty(sinks);
                            DEM = fillsinks(DEM);
                        else
                            DEM = fillsinks(DEM,sinks);
                        end
                        
                    case 'carve'
                        if verbose
                            disp([datestr(clock) ' -- Sink filling'])
                        end
                        if isempty(sinks);
                            DEMF = fillsinks(DEM);
                        else
                            DEMF = fillsinks(DEM,sinks);
                        end
                        
                        if isempty(weights);                            
                            D    = DEMF-DEM;
                        else
                            D    = -weights;
                        end
                        D = D.Z;
                        DEM  = DEMF;
                end

                % construct height graph
                if p.Results.internaldrainage
                    [Iobj,SILLSobj,IntBasin] = identifyflats(DEM);
                else
                    [Iobj,SILLSobj] = identifyflats(DEM);
                end
                
                I     = Iobj.Z;
                SILLS = SILLSobj.Z;
                clear Iobj Sillsobj
                
                % calculate sills for internal lake basins. These should be
                % located within the lake to force convergent flows
                if p.Results.internaldrainage
                    % There are various ways to get a lake center pixel
                    % 1. regional maxima of distance transform. Most
                    % elegant and good looking, but may produce several
                    % internal sills and requires much computation
                    % 2. randomly pick one pixel in each sink
                    STATS   = regionprops(IntBasin.Z,'PixelIdxList');
                    for r = 1:numel(STATS);
                        I(STATS(r).PixelIdxList(1)) = false;
                        SILLS(STATS(r).PixelIdxList(1)) = true;
                    end
                    clear IntBasin STATS
                end
                    
                
                if verbose
                    disp([datestr(clock) ' -- Flat sections identified'])
                end

                switch preprocess
                    case 'carve'
                        CarveMinVal = 0.1;
                        if ~isscalar(cweight);
                            D = (D + cweight);                            
                            D = linscale(D,0,100);
                        end
                 
                        CC = bwconncomp(I);                        
                        for r = 1:CC.NumObjects;
                            D(CC.PixelIdxList{r}) = (max(D(CC.PixelIdxList{r})) - D(CC.PixelIdxList{r})).^tweight + CarveMinVal;
                        end
                        
                        if ~isscalar(streams);
                        if islogical(streams);                        
                            D(streams) = CarveMinVal;
                        else
                            D = max(D - D.*min(streams,0.99999),CarveMinVal);
                        end
                        end
                        clear CC
                end
                

                dem = DEM.Z;
                % establish the connectivity between sills and flats
                [row,col] = find(SILLS);
                IXsill    = sub2ind(FD.size,row,col);
                rowadd = [-1 -1 0 1 1  1  0 -1];
                coladd = [ 0  1 1 1 0 -1 -1 -1];
                PreSillPixel = [];
                for r = 1:8;
                    rowp = row + rowadd(r);
                    colp = col + coladd(r);
                    
                    ValidRowColPair    = rowp>0 & colp>0 & rowp<=FD.size(1) & colp<=FD.size(2);
                    IXPreSill = sub2ind(FD.size,rowp(ValidRowColPair),colp(ValidRowColPair));
                    
                    PreSillPixel = [PreSillPixel;...
                        IXPreSill((dem(IXsill(ValidRowColPair)) == dem(IXPreSill)) & I(IXPreSill))];   %#ok<AGROW>
                    
                end
                
                clear row col IXsill rowadd coladd ValidRowColPair IXPreSill

                % costs to route over flats
                I = ~I;
                switch lower(preprocess)
                    case {'fill','none'}
                        D = bwdist(I,'euclidean');
                        mask = inf(FD.size,class(D));
                        mask(I) = 0;
                        D = (imreconstruct(D+1,mask) - D)*FD.cellsize;
                    case 'carve'

                end
                if verbose
                    disp([datestr(clock) ' -- Weights for graydist calculated'])
                end

                % distance transform
                D(I) = inf;
                D = graydist(double(D),double(PreSillPixel),'q') + 1;
                D(I) = -inf;

                if verbose
                    disp([datestr(clock) ' -- Auxiliary topography in flats calculated'])
                end
                
                
                if mexx
                    % mexed version: 
                    % identifies steepest downward neighbors in the DEM and 
                    % the distance grid obtained from graydist and performs
                    % a topological sort
                    clear I S demf                     
                    D(isinf(D)) = inf;
                    D(SILLS)    = 0;
                    clear SILLS
                    dem = double(dem);
                    % call to steepest neighbor. 
                    SE = steepestneighbor_mex(dem,D,FD.cellsize);
                    if verbose
                       disp([datestr(clock) ' -- Steepest neighbor identified'])
                    end
                    clear D
                    I  = isnan(dem);
                    clear dem
                    IX = find(SE==0 & ~I);
                    MV = nnz(I);
                    [FD.ix,FD.ixc] = tsort_mex(SE,uint32(IX),uint32(MV));
                    I  = FD.ix==0;
                    FD.ix(I) = [];
                    FD.ixc(I) = [];
                    
                else

                    % Find steepest neighbor
                    % sort pixels in D
                    % adapted from sortrows.m
                    clear PreSillPixel
                    [~,IXSortedFlats] = sort(D(:),'descend');
                    clear D

                    if verbose
                        disp([datestr(clock) ' -- Pixels sorted (1)'])
                    end

                    ndx = (uint32(1):uint32(nrc))';
                    ndx = ndx(IXSortedFlats);
                    clear IXSortedFlats

                    [~,FD.ix] = sort(dem(ndx),'descend');

                    if verbose
                        disp([datestr(clock) ' -- Pixels sorted (2)'])
                    end

                    FD.ix = uint32(FD.ix);
                    FD.ix = ndx(FD.ix);
                    clear ndx;



                    % a fast solution that has quite much memory overhead...
                    pp = zeros(FD.size,'uint32');
                    IX = (uint32(1):uint32(numel(dem)))';
                    pp(FD.ix) = IX;

                    % cardinal neighbors
                    IXC1 = imdilate(pp,[0 1 0; 1 1 1; 0 1 0]>0);
                    xxx1 = IXC1;
                    IX   = IXC1(FD.ix);
                    IXC1 = FD.ix(IX);
                    G1   = (dem(FD.ix)-dem(IXC1))/(FD.cellsize);
                    G1(FD.ix == IXC1) = -inf;

                    % diagonal neighbors
                    IXC2 = imdilate(pp,[1 0 1; 0 1 0; 1 0 1]>0);
                    xxx2 = IXC2;
                    IX   = IXC2(FD.ix);
                    IXC2 = FD.ix(IX);
                    G2   = (dem(FD.ix)-dem(IXC2))/(norm([FD.cellsize,FD.cellsize]));

                    % choose the steeper one

                    I    = G1<=G2 & xxx2(FD.ix)>xxx1(FD.ix);
                    FD.ixc = IXC1;
                    FD.ixc(I) = IXC2(I);

                    I = FD.ixc == FD.ix;
                    FD.ix(I) = [];
                    FD.ixc(I) = [];

                    % remove nans
                    I = isnan(dem);
                    FD.ixc(I(FD.ix)) = [];
                    FD.ix(I(FD.ix)) = [];
                end

                if verbose
                    disp([datestr(clock) ' -- Ordered topology established'])
                end
            end
        end
    end


    function FD = saveobj(FD)
        FD.fastindexing = false;
    end
    
    function FD = set.fastindexing(FD,val)
        
        validateattributes(val,{'numeric','logical'},{'scalar'})
        
        if val
            FD.fastindexing = true;
            FD.ixcix  = zeros(FD.size,'uint32');
            FD.ixcix(FD.ix) = uint32(1):uint32(numel(FD.ix));
        else
            FD.fastindexing = false;
            FD.ixcix  = [];
        end
    end

end
end


