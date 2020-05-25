function [CJ,varargout] = jctcon(D,varargin)
%JCTCON        compute junction connectivity
%
% Syntax
%
%     CJ = jctcon(D)
%     CJ = jctcon(D,maxdist)
%     [CJ,x,y] = jctcon(D)
%     
%     
% Description
%     
%     JCTCON computes the junction connectivity for junctions provided by
%     the linear indices in ixjct. The junction connectivity is defined to
%     be the sum of the ratios of the Euclidean distance, d , and the 
%     divide distance, d_d , of all divide edges, n, within a specified 
%     maximum divide distance, d_{d,max}. See Scherler and Schwanghart
%     (2020) for more details.
%     
%
% Input arguments
%
%     D       instance of DIVIDEobj
%     ixjct   linear index of junction(s)
%     maxdist maximum divide distance (default = 5000)
%     
% 
% Output arguments
%
%     CJ      Junction connectivity values
%     x,y     x,y coordinates of junctions 
%
% Examples
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','c');
%     ST = STREAMobj(FD,flowacc(FD)>1000);
%     D = DIVIDEobj(FD,ST);
%     D = divorder(D,'topo');
%     [CJ,x,y] = jctcon(D);
%     plot(D,'color',[.7 .7 .7])
%     hold on
%     scatter(x,y,50,CJ,'filled')
%     hc = colorbar;
%     hc.Label.String = 'Junction connectivity';
%
%
% See also: DIVIDEobj
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: April 2020


maxdist = 5e3;

if nargin==3
    maxdist = varargin{1};
end



CJ = nan(size(D.jct));

for i = 1 : length(D.jct) 
    
    D2 = resortdivide(D,D.jct(i),maxdist*2);
    D2 = divdist(D2);
    
    [x0,y0] = ind2coord(D,D.jct(i));
    [x,y] = ind2coord(D2,D2.IX);
    dist = hypot(x-x0,y-y0); % euclidean distance
    dd = D2.distance; % divide distance
    [dd,six] = sort(dd,'ascend');
    dist = dist(six);

    ix = dd<maxdist;
    d = dist(ix);
    dd = dd(ix);
    
    ix = d>0 & dd>0;
    r = d(ix)./dd(ix);
    CJ(i) = sum(r).*D.cellsize./maxdist;

end

if nargout>1
    [x,y] = ind2coord(D,D.jct);
    varargout{1} = x;
    varargout{2} = y;
end


end


function DOUT = resortdivide(DIN,ixep,distpenality)

DOUT = DIN;
DOUT.ep = ixep;

M = onl2struct(DIN.IX);
[M.flag] = deal(true);
[M.maxdist] = deal(nan);
segments = struct;
ia = max(ismember(vertcat(M.st),ixep),[],2);
[M(ia).dist] = deal(0);

maxdist = zeros(size(ixep));

ct = 0;
while ~isempty(ixep) 
    
    [ia,ib] = ismember(vertcat(M.st),ixep);
    [locr,locc] = find(ia);
    
    % flip segments, if needed
    ix = find(locc==2);
    for i = 1 : length(ix) 
        thisr = locr(ix(i));
        M(thisr).IX = [flipud(M(thisr).IX(1:end-1));NaN];
        M(thisr).st = fliplr(M(thisr).st);
    end
    
    % update distance
    for i = 1 : length(locr) 
        thisr = locr(i);
        thisc = locc(i);
        M(thisr).dist = maxdist(ib(thisr,thisc)) + (numel(M(thisr).IX)-2).*DIN.cellsize;
    end
    
    % loop over segments 
    ixep = [];
    maxdist = [];
    for i = 1 : length(locr)
        thisr = locr(i);
        if M(thisr).flag && M(thisr).dist<distpenality
            
            % Probably better to cut segments at max distance...
            % At present, only entire segments are considered.
            
            % transfer segment
            ct = ct+1;
            segments(ct).ix = M(thisr).IX;
            M(thisr).flag = false;
            ixep(ct) = M(thisr).st(2); % second, because flipped
            maxdist(ct) = M(thisr).dist;
        end
    end % for-endpoints loop
    
end % while-endpoints loop

DOUT.IX = vertcat(segments.ix);
DOUT.issorted = true;
DOUT.jct = [];
DOUT.jctedg = [];
DOUT.order = [];
DOUT.distance = [];


end


