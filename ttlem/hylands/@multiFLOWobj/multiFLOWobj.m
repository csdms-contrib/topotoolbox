classdef multiFLOWobj
    
    %multiFLOWobj create flow direction object
    %
    % Syntax
    %
    %     FD = multiFLOWobj(DEM)
    %     FD = multiFLOWobj(DEM,'pn','pv',...)
    %     FD = multiFLOWobj(DEM,'multi')
    %     FD = multiFLOWobj(DEM,'Dinf')
    %     FD = multiFLOWobj(M,'cellsize',cs,'size',siz)
    %
    % Description
    %
    %     multiFLOWobj creates a flow direction object which can be used to
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
    %     multiFLOWobj(DEM,'multi') derives multiple flow direction (MFD) and
    %     multiFLOWobj(DEM,'dinf') derives D infinity according to Tarboton's (1997)
    %     method (Eddins 2016). Note that this option precludes any other
    %     options to be set. Thus, this method needs preprocessing, e.g. using
    %     fillsinks(DEM).
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
    %     'preprocess' --  {'carve'}, 'fill', 'none'
    %            set DEM preprocessing that determines flow behavior in
    %            topographic depressions and flat areas
    %     'sinks' -- logical matrix same size as dem
    %            true values in sinks are treated as sinks in the digital
    %            elevation model and are not filled or carved, if the
    %            preprocessing option fill or carve are chosen.
    %     'internaldrainage' -- {false} or true
    %            set this parameter value to true if flow directions should be
    %            derived in the lowest, flat regions of internal drainage
    %            basins. By default, this parameter is set to false since this
    %            information is usually not required and flow paths will stop
    %            when entering flat, internally drained sections.
    %     'cweight' --  scalar {1}
    %            adjust cost matrix if preprocessing option 'carve' has been
    %            chosen.
    %     'verbose' --  {false},true
    %            verbose output in the command window to track computational
    %            progress. Particularly interesting when working with very
    %            large matrices.
    %     'mex' --  {'false'},true
    %            controls if the mex routines should be called. If true, the mex
    %            functions (see function compilemexfiles) must be available on
    %            the search path. Using mex functions can increase the speed at
    %            which an instance of multiFLOWobj is constructed.
    %
    % Applicable only, if calculated from flow direction matrix M
    %
    %     'cellsize' -- {1}
    %            cellsize in x and y direction (scalar double)
    %     'size' --  (1x2 array, e.g. size(dem))
    %            the size of the DEM of which the flow direction matrix was
    %            derived from.
    %     'refmat' -- referencing matrix as derived from GRIDobj property refmat
    %     'algorithm' -- {'dmperm'} or 'tsort' or 'toposort'
    %            determines the algorithm to perform a topological sorting
    %            of the vertices of the directed graph in matrix M (this option
    %            is left here mainly for evaluation purposes of different
    %            algorithms).
    %
    % Output
    %
    %     FD     flow direction object (multiFLOWobj)
    %
    % Example 1
    %
    %     % create and evaluate flow direction object for a DEM
    %     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
    %     FD = multiFLOWobj(DEM,'preprocess','carve');
    %     A  = flowacc(FD);
    %     imageschs(DEM,log(A));
    %
    % Example 2
    %
    %     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
    %     FD  = multiFLOWobj(DEM,'preprocess','carve','mex',true);
    %     DEM = imposemin(FD,DEM,0.0001);
    %     S   = STREAMobj(FD,'minarea',1000);
    %     DEMc = DEM;
    %     DEMc.Z(S.IXgrid) = DEMc.Z(S.IXgrid)-100;
    %     FD  = multiFLOWobj(DEMc,'multi');
    %     A   = flowacc(FD);
    %     imageschs(DEM,log(A));
    %
    %
    % See also: GRIDobj, STREAMobj
    %
    % References:
    %
    %     Tarboton, D. G. (1997). A new method for the determination of flow
    %     directions and upslope areas in grid digital elevation models.
    %     Water Resources Research, 33(2), 309-319.
    %
    %     Eddins, S. (2016). Upslope area function. Mathworks File Exchange,
    %     https://www.mathworks.com/matlabcentral/fileexchange/15818-upslope-area-functions
    %
    % Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
    % Date: 02. September, 2017
    
    
    % Update
    % 2015-06-16: added support for digraph toposort algorithm available with
    % R2015b
    % 2016-10-07: added support for multiple flow directions
    % 2016-11-14: added support for Dinf
    % 2017-09-02: default preprocessing option is carve
    
    properties(GetAccess = 'public', SetAccess = 'public')
        size      % size of instance of GRIDobj from which STREAMobj was derived
        type      % flow direction type (single, multi)
        ix        % [edge attribute] topologically sorted nodes (givers)
        ixc       % [edge attribute] topologically sorted nodes (receivers)
        fraction  % [edge attribute] fraction transfered between nodes
        cellsize  % cellsize of the grid (scalar)
        refmat    % 3-by-2 affine transformation matrix (see makerefmat)
        georef    % additional information on spatial referencing
    end
    
    properties(GetAccess = 'public', SetAccess = 'public')
        fastindexing = false; % set to true to initiate fast indexing
        ixcix     % indexing matrix for fast indexing
    end
    
    methods
        function FD = multiFLOWobj(DEM,varargin)
            
            
            %% Multiple flow direction
            FD.cellsize = DEM.cellsize;
            FD.refmat   = DEM.refmat;
            FD.georef   = DEM.georef;
            FD.size     = DEM.size;
            if numel(varargin)>1
                nb=varargin{2};
            end
            switch lower(varargin{1})
                case 'multi'
                    M = flowdir(DEM,'type','multi','nb',nb);
                    FD.type = 'multi';
                case 'dinf'
                    R = dem_flow(DEM.Z);
                    M = flow_matrix(DEM.Z,R,DEM.cellsize,DEM.cellsize);
                    M = -M';
                    
                    FD.type = 'Dinf';
                otherwise
                    error('unknown flow direction type')
            end
            
            try
                p = toposort(digraph(M));
                p = uint64(p);
            catch
                [p,~] = dmperm(speye(size(M))-M);
                p = uint64(p);
            end
            
            [FD.ix,FD.ixc,FD.fraction] = find(M(p,p));
            
            
            p = p(:);
            FD.ixc = p(FD.ixc);
            FD.ix  = p(FD.ix);
        end
    end
end