function dist = dist2node(DIN,ixep,varargin)
%DIST2NODE network distance to nodes
%
% Syntax
%
%     d = dist2node(DIN,IX)
%     d = dist2node(DIN,IX,dist)
%
% Description
%     
%     DIST2NODE calculates the distance along the divide network from
%     specific divide network nodes (endpoints or junctions). The input IX 
%     provides the linear indices to nodes in the divide network. If input 
%     dist is provided, only divide segments at a distance less than dist
%     will be considered.
%     
% Input arguments
%     
%     DIN    instance of divide network object (DIVIDEobj)
%     IX     linear indices to nodes of the input DIVIDEobj
%     dist   optional distance penality to restrict the extent of the 
%            output divide network
%     
% Output arguments
%
%     DOUT      instance of DIVIDEobj
%
% Examples
%     
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     FA = flowacc(FD);
%     ST = STREAMobj(FD,FA>500);
%     D = DIVIDEobj(FD,ST);
%     D = divorder(D,'topo');
%     ix = D.jct(1);
%     d = dist2node(D,ix);
%     plotc(D,d)
%  
% See also: DIVIDEobj, DIVIDEobj/sort
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: April 2020

distpenality = inf;
if nargin==3
    distpenality = varargin{1};
end


DOUT = DIN;
DOUT.ep = ixep;

II = (1:numel(DIN.IX))';

M = onl2struct(DIN.IX,'II',II);
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
        
        M(thisr).II = flipud(M(thisr).II);
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
DOUT = divdist(DOUT);


[~,locb] = ismember(DIN.IX,DOUT.IX);

dist = nan(size(locb));
ix = locb>0;
dist(ix) = DOUT.distance(locb(ix));


end

