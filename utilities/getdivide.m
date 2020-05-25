function OUT = getdivide(MS,IX,FD)
%GETDIVIDE drainage divide from drainage basin outline
%
% Syntax
%
%     MS2 = getdivide(MS,IX,FD)
%
% Description
%
%     Basin outlines obtained with the functions drainagebasins and 
%     GRIDobj2polygon are polygons that are defined by lines that start and
%     end at the drainage basin pixel with the lowest linear index of the
%     GRIDobj. This function takes these lines, and makes sure that they 
%     start and end at the stream outlets of the basins but don't cross 
%     them.
%
% Input arguments
%
%     MS      mapping structure (obtained from GRIDobj2polygon)
%     IX      linear indices to outlets (obtained from drainagebasins)
%     FD      instance of FLOWobj
%
%
% Output arguments
%
%     MS2     mapping structure
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','c');
%     FA = flowacc(FD);
%     ST = STREAMobj(FD,FA>1000);
%     S = streamorder(FD,FA>1000);
%     [DB,IX] = drainagebasins(FD,S,1);
%     MS = GRIDobj2polygon(DB);
%     MS = getdivide(MS,IX,FD);
%     plot(ST), hold on
%     plot(vertcat(MS.X),vertcat(MS.Y),'k-')
%     axis equal
%     
%
% See also: STREAMobj
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: September 2018


% The divide grid (DG) is shifted by half a cell size in x & y
cs = FD.cellsize;
hcs = cs/2;
fracs = hcs*1e-6;
DG = GRIDobj([]);
DG.size = FD.size+[1 1];
DG.refmat = FD.refmat+[0 0;0 0;-hcs,hcs];


% Identify internal drainage 
% get corner coordinates of outlets 
[tx,ty] = ind2coord(FD,IX);
fx = [tx-hcs,tx-hcs,tx+hcs,tx+hcs];
fy = [ty-hcs,ty+hcs,ty-hcs,ty+hcs];
fi = nan(size(fx));
fi(:) = coord2ind(DG,fx,fy);
allx = vertcat(MS.X);
ally = vertcat(MS.Y);
allix = coord2ind(DG,allx,ally);
isendo = not(any(ismember(fi,allix),2));
if sum(isendo)>0
    warning('Warning: found endorheic drainage. Results may be erroneous.')
end


%% Find nodes neighboring the outlets
[x1,y1] = ind2coord(FD,IX);
x2 = nan(size(IX)); y2 = x2;

% Downstream neighbors?
[isdo,dox] = ismember(IX,FD.ix);
if any(isdo) % tributary junctions with downstream neighbors
    [x2(isdo),y2(isdo)] = ind2coord(FD,FD.ixc(dox(isdo)));
end
if not(all(isdo)) % outlets without downstream neighbors
    % Create GRIDobj from FLOWobj with padded nan rim
    G = GRIDobj([]);
    pg = properties(G);
    pf = properties(FD);
    p  = intersect(pg,pf);
    for r = 1:numel(p)
        G.(p{r}) = FD.(p{r});
    end
    G.Z = nan(G.size,'single');
    %P = polygon2GRIDobj(G,MS); % Fails sometimes due to a bug
    % This is just a work-around to avoid a bug in polygon2GRIDobj
    P = pad(polygon2GRIDobj(pad(G,1),MS),-1);
    G.Z(P.Z) = 1;
    G = pad(G,1,nan);
    
    % Get 8 neighbors
    n = nnz(not(isdo));
    nbx = x1(not(isdo))+[-cs, 0, +cs, -cs, +cs, -cs, 0, +cs];
    nby = y1(not(isdo))+[+cs, +cs, +cs, 0, 0, -cs, -cs, -cs];
    nix = reshape(coord2ind(G,nbx,nby),n,8);
    NI = isnan(G.Z(nix));
    nb = nan(n,1);
    exo = sum(NI,2)>0;
    % select one neighbor
    for i = 1 : n
        if exo(i) % exorheic 
            tnb = nix(i,(NI(i,:)));
            p = randperm(length(tnb),1);
            nb(i) = tnb(p);
        end
    end
    [x2(not(isdo)),y2(not(isdo))] = ind2coord(G,nb);
end

% Determine neighborhood relationship and nodes
mx = mean([x1,x2],2);
my = mean([y1,y2],2);
di = coord2ind(DG,mx,my); % one diagonal node
hi = reshape(coord2ind(DG,[mx,mx],[y1+hcs,y1-hcs]),length(isdo),2); % two nodes
vi = reshape(coord2ind(DG,[x1+hcs,x1-hcs],[my,my]),length(isdo),2); % two nodes
dc = abs((x1-x2)./cs);
dr = abs((y1-y2)./cs);

MI = cell(size(IX));
isdiag = false(size(IX));
for i = 1 : length(IX) 
    if isnan(dc(i)) % endorheic
        MI{i} = nan;
    elseif dc(i)>fracs && dr(i)>fracs % diagonal neighbors
        MI{i} = di(i); 
        isdiag(i) = true;
    elseif dc(i)>fracs && dr(i)<fracs % horizontal neighbors at (mx,y1+/-hcs)
        MI{i} = hi(i,:);
    elseif dc(i)<fracs && dr(i)>fracs % vertical neighbors at (x1+/-hcs,my)
        MI{i} = vi(i,:);
    end
end


%% Loop over outlets and shift outlines 
OUT = MS;
for i = 1 : length(IX)
    ix = not(isnan(MS(i).X));
    tx = MS(i).X(ix); tx = tx(1:end-1);
    ty = MS(i).Y(ix); ty = ty(1:end-1);
    tix = coord2ind(DG,tx,ty);
    clear ix
    
    if not(isnan(MI{i})) % neighbor exists
        ix = sort(find(ismember(tix,MI{i})));
        if isdiag(i) % duplicate nodes for start and end points
            ix = ones(1,2).*ix;
        end
        if ix(1)~=1 || ix(2)~=nnz(tx) 
            tx = [tx(ix(2):end);tx(1:ix(1))];
            ty = [ty(ix(2):end);ty(1:ix(1))];
        end
        if isdiag(i) 
            tix = coord2ind(DG,tx,ty);
            aa = ismember(tix,MI{i});
            [tx,ty] = ind2coord(DG,tix(not(aa)));
        end
        %OUT(i).endo = false; %%%
        
    else % endorheic: do nothing
        %OUT(i).endo = true; %%%
        
    end
    
    OUT(i).X = [tx; nan];
    OUT(i).Y = [ty; nan];
    
end


end 