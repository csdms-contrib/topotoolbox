function d = distance(S,type)

%DISTANCE return node attribute list with distances along the stream network
%
% Syntax
%
%     d = distance(S,type)
%     d = distance(S,S2)
%     d = distance(S,IX)
%     d = distance(S,xy)
%
% Description
%
%     This function returns a node-attribute list with a distance value for  
%     each node in the stream network S. Since distance can be calculated
%     in different ways for a channelnetwork, the second input argument
%     lets you choose between a number of different distance metrics. If
%     the second input argument S2 is a STREAMobj, distances will be based
%     on the distances of S2. Note that S must be a subset of S2 then! 
%
%     If the second input argument is IX or xy then d will not be a
%     node-attribute list, but a vector with as many locations in IX or
%     xy.
%
% Input arguments
%
%     S       STREAMobj
%     type    'from_outlet'      distance in upstream direction
%             'min_from_ch'      shortest from channelhead
%             'max_from_ch'      longest from channelhead
%             'mean_from_ch'     mean distance from channelheads
%             'nr_of_ch'         number of channelheads (not a distance
%                                measure)
%             'node_to_node'     distance between each node and its
%                                downstream neighbor
%             'accumdownstream'  downstream accumulated distance
%             'from_trunk'       distance in upstream direction from trunk
%                                stream
%     S2      STREAMobj of which S is a subset.
%     IX      linear index into GRIDobj (points must be on the stream
%             network)
%     xy      n*2 matrix with coordinates. Points will be snapped to the
%             stream network (see STREAMobj/snap2stream). 
%
% Output arguments
%
%     d       node attribute list. If the second input argument is IX or xy
%             then d will be a vector with as many locations in IX or xy.
%
% Example 1: Planform stream network plot with colors showing distance from
%            channelheads
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM);
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S,1);
%     d = distance(S,'max_from_ch');
%     plotc(S,d/1000)
%     h = colorbar;
%     h.Label.String = 'Max. distance from channelhead [km]'; 
%
% Example 2: Plotting random points into a longitudinal profile plot
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM);
%     S = STREAMobj(FD,'minarea',1000);
%     IX = randlocs(S,100);
%     d = distance(S,IX);
%     plotdz(S,DEM)
%     hold on
%     plot(d,DEM.Z(IX),'s')
%     hold off
%
% Example 3: Plot a subset of a stream network using plotdz
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM);
%     S = STREAMobj(FD,'minarea',1000);
%     S2 = modify(S,'streamorder',2);
%     d = distance(S2,S);
%     plotdz(S,DEM)
%     hold on
%     plotdz(S2,DEM,'LineWidth',2,'distance',d);
%     hold off
%
% See also: STREAMobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 9. August, 2018

if nargin == 1
    d = S.distance;
    return
end

if ischar(type)
    validtypes = {...
        'from_outlet', ... % distance in upstream direction
        'min_from_ch',... % shortest from channelhead
        'max_from_ch',... % longest from channelhead
        'mean_from_ch',... % mean distance from channelheads
        'nr_of_ch'...
        'node_to_node'...
        'accumdownstream' ...
        'from_trunk' ...
        };
    
    type = validatestring(type,validtypes,'distance','type',2);
    
    % handle the simple case which is implement as dynamic property of S
    switch type
        case 'from_outlet'
            d = S.distance;
            return;
    end
    
    % preallocate array
    d = nan(size(S.x));
    
    % distance between two nodes
    dedge = sqrt((S.x(S.ixc)-S.x(S.ix)).^2 + (S.y(S.ixc)-S.y(S.ix)).^2);
    
    switch type
        case 'node_to_node'
            d(:) = 0;
            d(S.ix) = dedge;
            return
    end
    
    % identify channel heads
    M  = sparse(double(S.ix),double(S.ixc),true,numel(d),numel(d));
    CH = (sum(M,1)'==0) & (sum(M,2) ~= 0);
    
    % set distance at channel heads to zero
    d(CH) = 0;
    switch type
        case 'min_from_ch'
            for r = 1:numel(S.ix)
                d(S.ixc(r)) = min(d(S.ix(r)) + dedge(r),d(S.ixc(r)));
            end
        case 'max_from_ch'
            for r = 1:numel(S.ix)
                d(S.ixc(r)) = max(d(S.ix(r)) + dedge(r),d(S.ixc(r)));
            end
        case 'nr_of_ch'
            nrc = numel(d);
            d = full((speye(nrc)-M')\double(CH));
        case 'mean_from_ch'
            nrc = numel(d);
            nr = (speye(nrc)-M')\double(CH);
            
            d  = accumarray(S.ixc,dedge,[nrc,1],@mean,0).*nr;
            
            d  = (speye(nrc)-M')\d;
            d  = full(d./nr);
        case 'accumdownstream'
            d(:) = 0;
            for r = 1:numel(S.ix)
                d(S.ixc(r)) = d(S.ix(r))+dedge(r) + d(S.ixc(r));
            end
        case 'from_trunk'
            d(:) = 0;
            St   = trunk(S);
            I    = ~ismember(S.IXgrid,St.IXgrid);
            for r = numel(S.ix):-1:1
                if I(S.ix(r))
                    d(S.ix(r)) = d(S.ixc(r)) + dedge(r);
                end
            end
                    
    end
    
elseif isa(type,'STREAMobj')
    S2 = type;
    d2 = S2.distance;
    [I,inS2] = ismember(S.IXgrid,S2.IXgrid);
    if ~all(I)
        error('S is not a subset of S2');
    end
    d  = d2(inS2);
    
else
    if isscalar(type) || iscolumn(type)
        IX = type;
    else 
        [~,~,IX] = snap2stream(S,type(:,1),type(:,2));
    end
    
    [I,ix] = ismember(IX,S.IXgrid);
    if ~all(I)
        error('Some of the query points are not on the stream network')
    end
    
    d = S.distance;
    d = d(ix);
    
end
    
