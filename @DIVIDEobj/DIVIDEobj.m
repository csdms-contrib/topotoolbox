classdef DIVIDEobj
%DIVIDEobj Create divide object (DIVIDEobj)
%
% Syntax
%
%     D = DIVIDEobj(FD,ST) 
%     D = DIVIDEobj(FD,ST,pn,pv)
%
% Description
%
%     An instance of divide object encapsulates information on the geometry
%     and connectivity of a divide network, based on the flow direction of
%     a digital elevation model and a STREAMobj or a logical raster that 
%     indicates the position of streams. Divides are defined as the lines
%     surrounding drainage basins and thus they are positioned in between
%     pixels of the digital elevation model. The drainage basins used to
%     define divide objects are based on tributrary junctions that are
%     derived from changes in stream orders.
%     DIVIDEobj provides access to various methods that investigate 
%     properties of a divide network and associated data.
%
% Input arguments
%
%     FD     instance of flow direction object (FLOWobj)
%     ST     instance of stream object (STREAMobj) or logical grid 
%            (GRIDobj) that indicates stream locations (e.g. obtained from 
%            thresholding the flow accumulation raster)
%
% Parameter name/value pairs
%
%     type         string that defines the stream ordering scheme 
%                  ('strahler', 'shreve', or 'topo' {default}) used to
%                  locate tributary junctions at breaks in stream orders
%     outlets      toggle (default=true) to indicate whether outlets, i.e.,
%                  the downstream ends of streams, shall be used to derive 
%                  divides
%     verbose      toggle for displaying function execution progress in the
%                  command window
%
% Output arguments
%
%     D      instance of DIVIDEobj
%
% Examples
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','c');
%     ST = STREAMobj(FD,flowacc(FD)>1000);
%     D = DIVIDEobj(FD,ST);
%     plot(D)
%  
% See also: getdivide, FLOWobj/drainagebasins, FLOWobj/streamorder
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: August 2020

    
    properties 
        size      % size of hypothetical GRIDobj that would place the 
                  % divide nodes into cells
        cellsize  % cellsize
        refmat    % 3-by-2 affine transformation matrix (see makerefmat)
        IX        % nan-separated linear indices of divide nodes into 
                  % instance of GRIDobj
                  
        order     % divide order
        distance  % maximum directed distance from divide endpoint        
        ep        % linear indices of divide network endpoints 
        jct       % linear indices of divide network junctions
        jctedg    % number of edges per junction (divides)
        %endo      % linear indices of endorheic nodes
        issorted  % flag to indicate if divide network is sorted
        ordertype % name of ordering scheme
    end
    
    
    methods 
        
        function [D,varargout] = DIVIDEobj(FD,ST,varargin) 
            %DIVIDEobj Construct an instance of this class
            
            % Parse inputs
            ordschemes = {'strahler','shreve','topo'};
            p = inputParser;
            p.FunctionName = 'divides';
            p.KeepUnmatched = true;
            addRequired(p,'FD',@(x) isa(x,'FLOWobj'));
            addRequired(p,'ST',@(x) isa(x,'GRIDobj')|isa(x,'STREAMobj'));
            addParamValue(p,'type','topo',@(x) ismember(lower(x),ordschemes))
            addParamValue(p,'outlets',true,@(x) islogical(x))
            addParamValue(p,'network',true,@(x) islogical(x))
            addParamValue(p,'verbose',false,@(x) islogical(x))
            parse(p,FD,ST,varargin{1:end});
            
            type    = p.Results.type;
            outlets = p.Results.outlets;
            verbose = p.Results.verbose;
            
            % Prepare
            % the divide grid is shifted by half a cell size in x & y
            hcs = FD.cellsize/2;
            D.cellsize = FD.cellsize;
            D.size = FD.size+[1 1];
            D.refmat = FD.refmat+[0 0;0 0;-hcs,hcs];
            
            % Create streamorder grid
            if isa(ST,'STREAMobj')
                STG = STREAMobj2GRIDobj(ST);
            elseif isa(ST,'GRIDobj')
                STG = ST;
            else
                error('TT2: See function help for instructions.')
            end
            sotype = lower(type);
            if strcmp(sotype,'topo')
                sotype = 'shreve';
            end
            S = streamorder(FD,STG,sotype);
            sz = S.Z(S.Z>0);
            sox = unique(sz);
            n = length(sox);
            
            % Divide identification loop
            G = cell(n+double(outlets),1); % divides
            O = cell(n+double(outlets),1); % outlets
            if verbose
                fprintf(1,'Divide identification:\n');
            end
            for k = 1 : n
                if verbose
                    fprintf(1,'%d / %d  \n',k,n);
                end
                [DB,IX] = drainagebasins(FD,S,sox(k));
                % Record divide coordinates
                if sum(IX)>0
                    % Get basin outlines
                    MS = GRIDobj2polygon(DB);
                    MS = getdivide(MS,IX,FD);
                    G{k} = MS;
                    O{k} = double(IX);
                end
            end
            if (outlets) % Include outlets
                IX = streampoi(ST,'outlets','ix');
                DB = drainagebasins(FD,IX);
                [MS,~,~] = GRIDobj2polygon(DB);
                MS = getdivide(MS,IX,FD);
                G{k+1} = MS;
                O{k+1} = IX;
            end
            
            
            % Make structure (M) with all divides
            M = struct;
            ct = 0;
            for k = 1 : length(G)
                MS = G{k};
                for j = 1 : length(MS)
                    x = MS(j).X;
                    y = MS(j).Y;
                    if length(x)>2 % At least two nodes and one NaN
                        ct = ct+1;
                        M(ct).IX = coord2ind(D,x,y);
                        M(ct).ep = M(ct).IX([1 end-1]);
                        %M(ct).endo = MS(j).endo.*ones(size(M(ct).IX));
                    end
                end
            end
            I = vertcat(M.IX); % nodes
            %D.ep = vertcat(M.ep); % all divide endpoints
            
            % remember endorheic nodes
            %ex = logical(vertcat(M.endo));
            %endo = unique(I(ex));
            
            % Remove redundant edges
            T = [I(1:end-1),I(2:end)]; % edges (connections between nodes)
            % Rows with one NaN in either column become NaN
            T(sum(isnan(T),2)==1,:) = NaN;
            % Set redundant edges to NaN
            sT = sort(T,2); % lower index left
            [~,IX,~] = unique(sT,'rows','stable');
            T0 = nan(size(T));
            T0(IX,:) = T(IX,:);
            % Remove redundant NaN between segments
            IX = isnan(T0(:,1));
            JX = [IX(2:end);true];
            IJ = logical(min([IX JX],[],2));
            T0 = [T0(not(IJ),:);NaN NaN];
            
            % Assemble structure
            M = struct;
            ix = [0;find(isnan(T0(:,1)))];
            for i = 1 : length(ix)-1 
                M(i).IX = [T0(ix(i)+1,1);T0(ix(i)+1:ix(i+1),2)];
                % % identify endorheic nodes
                % M(i).isendo = ismember(M(i).IX,endo);
            end
            D.IX = vertcat(M.IX);
            %D.endo = endo;
            D.issorted = false;
            
            % Get junctions and endpoints
            if p.Results.network
                D = divnet(D,FD);
                D = sort(D);
                D = divdist(D);
            end
            
            if nargout>1
                varargout{1} = O;
            end
            
        end
        
        
        function DOUT = divnet(DIN,FD) 
            %DIVNET   Compute divide network 
            % 
            % Syntax
            %
            %     D2 = divnet(D,FD)
            %
            % Description
            %
            %     DIVNET finds the endpoints and junctions that define the
            %     divide network in a divide object. Endpoints are the
            %     outer limits of a divide network, equivalent to channel
            %     heads in a drainage network. Junctions are points where
            %     two (or more) divides meet. They are equivalent to
            %     confluences in a drainage network. Both endpoints and
            %     junctions will be stored in the form of linear indices in
            %     the divide object.
            %
            % Input
            %
            %     D         instance of class DIVIDEobj
            %     FD        instance of class FLOWobj
            %
            % Output
            %
            %     D2         instance of class DIVIDEobj
            %
            % Example
            %
            %     D = divnet(D,FD);
            %
            % See also: DIVIDEobj, DIVIDEobj/sort
            %
            % Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
            % Date: Nov 2018
            
            %  input properties: D.IX
            % output properties: D.ep, D.jct, D.jctedg
            
            DOUT = DIN;
            
            % Grid with nodes indicating diagonal flow
            [x,y] = refmat2XY(FD.refmat,FD.size);
            cs = FD.cellsize;
            hcs = cs/2;
            xn = [x-hcs,x(end)+hcs];
            yn = [y+hcs;y(end)-hcs];
            A = zeros(FD.size+1);
            FX = GRIDobj(xn,yn,A);
            NEDGE = FX;
            NST = FX;
            [x,y] = getcoordinates(FX);
            [XD,YD] = meshgrid(x,y);
            
            % Indices of flow edge midpoints 
            [x1,y1] = ind2coord(FD,FD.ix);
            [x2,y2] = ind2coord(FD,FD.ixc);
            dx = x2-x1;
            dy = y2-y1;
            mx = (x1+x2)./2;
            my = (y1+y2)./2;
            [mix,res] = coord2ind(DIN,mx,my);
            FX.Z(mix(res<0.5*hcs & abs(dx+dy)<cs)) = 1; % NW-SE / SE-NW
            FX.Z(mix(res<0.5*hcs & abs(dx+dy)>cs)) = 2; % NE-SW / SW-NE
            
            % Unique nodes and number of edges
            T0 = [DIN.IX(1:end-1),DIN.IX(2:end)];
            T0(sum(isnan(T0),2)==1,:) = NaN;
            ix = not(isnan(T0(:,1)));
            T1 = T0(ix,:);
            [e,~,f] = unique([T0(ix,1);T0(ix,2)],'rows','stable');
            CTS = [e accumarray(f,1)];
            NEDGE.Z(CTS(:,1)) = CTS(:,2);
            
            % Unique segment termini and number of appearances
            M = onl2struct(DIN.IX);
            st = [M.st];
            [e,~,f] = unique(st(:),'stable');
            STIX = [e accumarray(f,1)];
            NST.Z(STIX(:,1)) = STIX(:,2);
            
            % Junctions
            DOUT.jct = find(NEDGE.Z>2 & FX.Z==0);
            DOUT.jctedg = NEDGE.Z(DOUT.jct);
            
%             % Endpoints
%             ixep = find((NEDGE.Z==1 | NEDGE.Z==3) & NST.Z==1);
                
            % Dead segments
            ixdseg = find((NEDGE.Z==2 & NST.Z==2 & FX.Z==0) | ...
                         (NEDGE.Z==4 & NST.Z==2 & FX.Z>0) | ...
                         (NEDGE.Z==4 & NST.Z==4 & FX.Z>0) | ...
                         (NEDGE.Z==3 & NST.Z==3 & FX.Z>0));
                     
            % Distinguish EP and DST at NST==2,NEDGE==2,ST>0
            ix = find(NEDGE.Z==2 & NST.Z==2 & FX.Z>0);
            L = ix;
            for i = 1 : length(ix)
                [r,c] = find(T1==ix(i));
                x1 = XD(T1(r(1),3-c(1))); % 3-c=2 if c=1 and 3-c=1 if c=2
                y1 = YD(T1(r(1),3-c(1)));
                x2 = XD(T1(r(2),3-c(2)));
                y2 = YD(T1(r(2),3-c(2)));
                dx = x2-x1;
                dy = y2-y1;
                L(i) = (abs(dx+dy)>cs)+1;
            end
            newix = not(logical(abs(L-FX.Z(ix))));
            %ixep = [ixep; ix(not(newix))];
            ixdseg = [ixdseg; ix(newix)];
            
            % Merge dead segment termini 
            M = onl2struct(DIN.IX);
            for i = 1 : length(M)
                ix = ismember(M(i).st,ixdseg);
                M(i).dead_st = ix;
            end
            allst = vertcat(M.st);
            dst = allst(vertcat(M.dead_st));
            udst = unique(dst(:));
            [M.tag] = deal(true);
            nn = length(udst);
            for i = 1:length(udst) 
                this_dst = udst(i);
%                 if verbose
%                     fprintf(1,'%1.0d / %1.0d\n',i,nn)
%                 end
                [r,c] = find(allst==this_dst);
                if length(r)==2 % 2 segment termini 
                    ix1 = M(r(1)).IX(1:end-1);
                    ix2 = M(r(2)).IX(1:end-1);
                    if c(1) == 1
                        ix1 = flip(ix1);
                    end
                    if c(2) == 2
                        ix2 = flip(ix2);
                    end
                    % Merge segments
                    M(r(1)).IX = [ix1;ix2(2:end);NaN];
                    M(r(1)).st = M(r(1)).IX([1,end-1])';
                    M(r(2)).tag = false;
                    
                elseif length(r)>2 % 3 or 4 segment termini 
                    k = FX.Z(this_dst);
                    % Get adjacent nodes
                    [tr,tc] = ind2sub(DIN.size,this_dst);
                    switch k
                        case 1 % 1 = NW-SE
                            p1 = [tr-1,tc;tr,tc+1];
                            p2 = [tr,tc-1;tr+1,tc];
                        case 2 % 2 = NE-SW
                            p1 = [tr-1,tc;tr,tc-1];
                            p2 = [tr,tc+1;tr+1,tc];
                    end
                    p{1} = sub2ind(DIN.size,p1(:,1),p1(:,2));
                    p{2} = sub2ind(DIN.size,p2(:,1),p2(:,2));
                    % Find and merge segments (once or twice)
                    for g = 1 : 2 
                        m = nan(1,2); ct = 0;
                        for j = 1 : length(r)
                            ax = ismember(p{g},M(r(j)).IX);
                            if any(ax)
                                ax = find(ax,1);
                                ct = ct+1;
                                m(ct) = j;
                                p{g}(ax) = [];
                            end
                        end
                        if isempty(p{g}) %not(any(isnan(m)))
                            ix1 = M(r(m(1))).IX(1:end-1);
                            ix2 = M(r(m(2))).IX(1:end-1);
                            if c(m(1)) == 1
                                ix1 = flip(ix1);
                            end
                            if c(m(2)) == 2
                                ix2 = flip(ix2);
                            end
                            % Merge segments
                            M(r(m(1))).IX = [ix1;ix2(2:end);NaN];
                            M(r(m(1))).st = M(r(m(1))).IX([1,end-1])';
                            M(r(m(2))).tag = false;
                            if length(r)==3
                                break;
                            end
                        end
                    end
                    
                end
                % Remove redundant segments 
                tags = [M.tag];
                M = M(tags);
                allst = vertcat(M.st);
            end
            M = onl2struct(vertcat(M.IX)); % Update structure

            % Insert breaks at junctions
            for i = 1:length(M) 
                ix = M(i).IX(1:end-1); 
                tix = ix;
                tix([1 end]) = NaN; % omit segment termini
                iy = find(ismember(tix,DOUT.jct));
                if nnz(iy)>0
                    iy = [1;iy(:);length(ix)];
                    NM = struct;
                    for k = 1 : length(iy)-1
                        NM(k).IX = [ix(iy(k):iy(k+1));NaN];
                    end
                    M(i).IX = vertcat(NM.IX);
                end
            end
            M = onl2struct(vertcat(M.IX)); % Update structure
            
            DOUT.IX = vertcat(M.IX);
            st = [M.st];
            DOUT.ep = setdiff(st(:),DOUT.jct);
            
        end
        
        
    end
end

