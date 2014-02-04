function OUT = tidyswath(SW,varargin)
% remove overlapping points from SWATHobj
%
% Syntax
%
%     OUT = tidyswath(SW)
%     OUT = tidyswath(SW,FD)
%     OUT = tidyswath(SW,FD,'both')
%
% Description
%
%     TIDYSWATH(SW) uses the function 'lineSegmentIntersect' by U. Murat 
%     Erdem, freely available on Matlab Central, to remove overlapping 
%     points in a SWATHobj by setting the corresponding Z values to NaN.
%
%     TIDYSWATH(SW,FD) assumes the SWATHobj SW was generated with a 
%     STREAMobj derived from the FLOWobj FD and sets data points outside 
%     the corresponding watershed to NaN.
%
%     TIDYSWATH(SW,FD,'both') does both of the above.
%
% Input arguments
%
%     SW     instance of SWATHobj
%     FD     instance of FLOWobj
%
% Output arguments
%
%     OUT    instance of SWATHobj
%
% Examples
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     FA  = flowacc(FD);
%     S = STREAMobj(FD,FA>1000);
%     S = klargestconncomps(S,1);
%     S = removeshortstreams(S,1e3);
%     SW = SWATHobj(DEM,S,'smooth',300,'plot',false);
%     SW = tidyswath(SW,FD,'both');
%     plot(SW)
%
% Author: Dirk Scherler (scherler[at]caltech.edu)
% Date: June, 2013



narginchk(1,3)

if nargin>=2
    FD = varargin{1};
    if ~isa(FD,'FLOWobj')
        error('Second input needs to be of class FLOWobj')
    end
end
if nargin==3
    if ~strcmp('both',varargin{2})
        error('Unknown third input to function call')
    end
end

if ~isa(SW,'SWATHobj')
    error('First input needs to be of class SWATHobj')
end


OUT = SW;

if nargin==1 || nargin==3
    
    % remove duplicates from mapping
    if ~exist('lineSegmentIntersect','file')
        error('Need function lineSegmentIntersect. Download from Matlab Central');
    else
        % loop of items in SWATHobj
        for k = 1 : length(SW.xy0)
            [ny,~] = size(SW.X{k});
            if mod(ny,2)==0; ctr = ny/2; xtr = 1;
            else ctr = (ny+1)/2; xtr = 0; end

            % do each side of the SWATHobj individually
            for sides = 1 : 2
                if sides==1; edge = 1;
                else edge = ny; ctr = ctr+1+xtr; end
                if edge < ctr; step = -1; else step = 1; end
                id = ctr:step:edge;
                d = abs(SW.disty{k}(id));
                maxd = max(d);

                % the points near the center of the swath are defined as the starting points
                x1 = SW.X{k}(ctr,:)';
                y1 = SW.Y{k}(ctr,:)';
                % the points on the edge of the swath are defined as the ending points
                x2 = SW.X{k}(edge,:)';
                y2 = SW.Y{k}(edge,:)';
                XY1 = [x1,y1,x2,y2];
                XY1(isnan(XY1(:,1)),:) = [];
                XY2 = XY1;
                LSI = lineSegmentIntersect(XY1,XY2);
                U = triu(LSI.intAdjacencyMatrix);

                % check each intersection
                [nru,ncu] = size(U);
                IX = find(U);
                for i = 1 : length(IX)
                    d1 = abs( LSI.intNormalizedDistance1To2(IX(i)) );
                    d2 = abs( LSI.intNormalizedDistance2To1(IX(i)) );
                    [p1,p2] = ind2sub([nru,ncu],IX(i));
                    % determine which line to cut short
                    if d1 <= d2;
                        thisp = p2;
                        thisd = d2*maxd;
                    else
                        thisp = p1;
                        thisd = d1*maxd;
                    end
                    % exclude points along line beyond intersection
                    ix = d>thisd;
                    OUT.Z{k}(id(ix),thisp) = nan;
                    %OUT.X{k}(id(ix),thisp) = nan;
                    %OUT.Y{k}(id(ix),thisp) = nan;
                end
            end
        end
    end

end

if exist('FD','var')
    
    FA = flowacc(FD);
    
    % collect outlets
    outl = zeros(size(SW.Z));
    outl_fa = zeros(size(SW.Z));
    for i = 1:length(SW.Z)
        x = SW.xy0{i}(:,1);
        y = SW.xy0{i}(:,2);
        ix = coord2ind(FA,x,y);
        fa = FA.Z(ix);
        [~,ixfa] = sort(fa,'descend');
        ix = ix(ixfa);
        outl(i) = ix(2);
        outl_fa(i) = fa(2);
    end
    
    [~,oix] = sort(outl_fa,'ascend');
    outl_s = outl(oix);

    ML = FLOWobj2GRIDobj(FD);
    ML.Z(:) = false;
    for i = 1 : length(outl_s)
        this_i = oix(i);
        L = drainagebasins(FD,outl_s(i));

        ix_sw = 1:numel(OUT.X{this_i});
        x = OUT.X{this_i}(ix_sw);
        y = OUT.Y{this_i}(ix_sw);

        % set out-of-range subscripts to nan
        xlims = FD.georef.SpatialRef.XLimWorld;
        ylims = FD.georef.SpatialRef.YLimWorld;
        ix = find(x>=xlims(1) & x<=xlims(2) & ...
            y>=ylims(1) & y<=ylims(2));
        ix_sw = ix_sw(ix);
        x = x(ix);
        y = y(ix);

        ix_grid = coord2ind(FD,x,y);
        ix_sw   = ix_sw(~isnan(ix_grid));
        ix_grid = ix_grid(~isnan(ix_grid));

        % set out-of-watershed to nan
        nix = ~L.Z(ix_grid);
        OUT.Z{this_i}(ix_sw(nix)) = nan;
        %OUT.Y{this_i}(ix_sw(nix)) = nan;
        %OUT.X{this_i}(ix_sw(nix)) = nan;

        nix = logical(ML.Z(ix_grid));
        OUT.Z{this_i}(ix_sw(nix)) = nan;
        %OUT.Y{this_i}(ix_sw(nix)) = nan;
        %OUT.X{this_i}(ix_sw(nix)) = nan;
        ML = ML|L;

    end

    
end

